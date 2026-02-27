if (!requireNamespace("reticulate", quietly = TRUE)) install.packages("reticulate")
library(reticulate)

ensure_pydeps= function() {
  invisible(reticulate::py_available(initialize = TRUE))
  for (pkg in c("numpy","pandas")) {
    if (!reticulate::py_module_available(pkg)) reticulate::py_install(pkg, pip = TRUE)
  }
  if (!reticulate::py_module_available("cv2")) reticulate::py_install("opencv-python-headless", pip = TRUE)
}

# --- simple BGR helper (optional; seeds usually keep one color) ---
.resolve_color_bgr= function(col) {
  if (is.numeric(col) && length(col) == 3) return(as.integer(col))
  if (!is.character(col) || length(col) != 1) stop("outline_color must be name or c(B,G,R)")
  col= tolower(col)
  pal= list(
    green = c(0L,255L,0L), red = c(0L,0L,255L), blue = c(255L,0L,0L),
    yellow = c(0L,255L,255L), white = c(255L,255L,255L),
    gray = c(180L,180L,180L), black = c(0L,0L,0L)
  )
  if (!col %in% names(pal)) stop("Unknown outline_color: ", col)
  pal[[col]]
}
#' Grain Morphometric Analysis from Digital Images
#'
#' This function performs high-throughput grain phenotyping by extracting physical traits
#' (area, perimeter, length, width) from images using OpenCV-based segmentation.
#'
#' @param input_folder Character. Path to the directory containing raw grain images (JPG, PNG).
#' @param output_folder Character. Path where processed images and the CSV result will be saved.
#' @param image_real_cm Numeric vector. Physical dimensions of the captured image area, e.g., c(width, height) in cm. Default is c(20, 20).
#' @param lower_hsv Integer vector. Lower bound for HSV color segmentation, e.g., c(H, S, V).
#' @param upper_hsv Integer vector. Upper bound for HSV color segmentation, e.g., c(H, S, V).
#' @param extra_hsv List of lists. Additional HSV ranges to include (useful for dark hilums or patterned seeds).
#' Each sub-list should contain 'lower' and 'upper' vectors.
#' @param margin Integer. Number of pixels to exclude from the image borders.
#' @param k_open Integer vector. Kernel size for morphological 'Opening' to remove background noise (dust).
#' @param k_close Integer vector. Kernel size for morphological 'Closing' to bridge cracks or gaps within seeds.
#' @param open_iter Integer. Number of iterations for the opening operation.
#' @param close_iter Integer. Number of iterations for the closing operation.
#' @param fill_holes Logical. If TRUE, force-fills internal voids/holes within detected seeds.
#' @param min_component_area_px Integer. Minimum pixel area for a candidate object to be considered.
#' @param object_min_area_cm2 Numeric. Minimum physical area (cm^2) to filter out non-grain debris.
#' @param rel_min_frac_of_largest Numeric. Relative filter (0 to 1); ignores objects smaller than this fraction of the largest detected seed.
#' @param max_keep Integer. Maximum number of grain objects to extract per image.
#' @param outline_color Character. Color of the detection boundary in output images (e.g., "green", "red").
#' @param outline_thickness Integer. Thickness of the contour lines in the processed output images.
#' @param ignore_processed_images Logical. If TRUE, skips files that already contain "_processed" in their filename.
#'
#' @returns A data frame containing morphometric data (Area, Perimeter, Width, Length) for each detected grain.
#' @export
#'
#' @examples
#' \dontrun{
#' # to install phenokio package
#' if(!require(remotes)) install.packages("remotes")
#' if (!requireNamespace("phenokio", quietly= TRUE)) {
#'  remotes::install_github("agronomy4future/phenokio", force= TRUE)
#' }
#' library(remotes)
#' library(phenokio)
#'
#' # Basic usage for soybean analysis
#' phenokio(
#'   input_folder= "path/to/images",
#'   output_folder= "path/to/results",
#'   image_real_cm= c(23, 23),
#'   lower_hsv= c(10, 80, 70),
#'   upper_hsv= c(35, 255, 255),
#'   fill_holes= TRUE,
#'   outline_color= "red"
#' )
#' * Github: https://github.com/agronomy4future/phenokio
#' # All Rights Reserved © J.K Kim (kimjk@agronomy4future.com)
#' }
phenokio= function(
    input_folder,
    output_folder,
    image_real_cm = c(20, 20),

    # --- segmentation (HSV) ---
    lower_hsv = c(0L, 0L, 0L),
    upper_hsv = c(180L, 255L, 255L),
    extra_hsv = NULL,              # keep for rare cases (two-tone seeds)

    # --- border exclusion ---
    margin = 0L,

    # --- morphology (most useful knobs for seeds) ---
    k_open  = c(5L, 5L),
    k_close = c(7L, 7L),
    open_iter  = 1L,
    close_iter = 2L,
    fill_holes = FALSE,            # seeds: often FALSE; set TRUE if dark pits split

    # --- filtering ---
    min_component_area_px = 200L,
    object_min_area_cm2 = 0.01,
    rel_min_frac_of_largest = 0.00, # 0 disables relative filter; use 0.05~0.15 if noisy
    max_keep = 2000L,

    # --- output overlay style (minimal) ---
    outline_color = "green",
    outline_thickness = 2L,

    ignore_processed_images = TRUE
) {
  ensure_pydeps()

  if (length(image_real_cm) == 1) image_real_cm <- rep(image_real_cm, 2)
  dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

  outline_color_bgr= .resolve_color_bgr(outline_color)

  py_code= "
import cv2
import numpy as np
import os, pandas as pd
from pathlib import Path
import re

def _as_uint8_triplet(x):
    arr = np.array(x, dtype=np.int64).reshape(-1)
    if arr.size != 3:
        raise ValueError('HSV must have 3 elements (H,S,V).')
    return np.clip(arr, 0, 255).astype(np.uint8)

def _combine_hsv_masks(hsv_img, lower_hsv, upper_hsv, extra_hsv):
    base = cv2.inRange(hsv_img, _as_uint8_triplet(lower_hsv), _as_uint8_triplet(upper_hsv))
    if extra_hsv is None:
        return base
    masks = [base]
    try:
        iterable = list(extra_hsv)
    except Exception:
        iterable = []
    for item in iterable:
        lo = hi = None
        try:
            lo = item.get('lower', None); hi = item.get('upper', None)
        except Exception:
            try:
                lo = item[0]; hi = item[1]
            except Exception:
                pass
        if lo is None or hi is None:
            continue
        m = cv2.inRange(hsv_img, _as_uint8_triplet(lo), _as_uint8_triplet(hi))
        masks.append(m)

    out = masks[0]
    for m in masks[1:]:
        out = cv2.bitwise_or(out, m)
    return out

def _apply_margin(mask, margin):
    m = int(margin) if margin is not None else 0
    if m <= 0:
        return mask
    out = mask.copy()
    h, w = out.shape[:2]
    out[:m, :] = 0
    out[h-m:h, :] = 0
    out[:, :m] = 0
    out[:, w-m:w] = 0
    return out

def _morph_open_close(mask, k_open, k_close, open_iter, close_iter):
    out = mask
    ko = tuple(int(x) for x in k_open)
    kc = tuple(int(x) for x in k_close)
    if ko[0] > 0 and ko[1] > 0 and int(open_iter) > 0:
        ker_o = np.ones(ko, np.uint8)
        out = cv2.morphologyEx(out, cv2.MORPH_OPEN, ker_o, iterations=int(open_iter))
    if kc[0] > 0 and kc[1] > 0 and int(close_iter) > 0:
        ker_c = np.ones(kc, np.uint8)
        out = cv2.morphologyEx(out, cv2.MORPH_CLOSE, ker_c, iterations=int(close_iter))
    return out

def _fill_holes(mask):
    inv = cv2.bitwise_not(mask)
    h, w = mask.shape[:2]
    flood = inv.copy()
    ffmask = np.zeros((h+2, w+2), dtype=np.uint8)
    cv2.floodFill(flood, ffmask, (0,0), 255)
    flood_inv = cv2.bitwise_not(flood)
    return cv2.bitwise_or(mask, flood_inv)

def _is_processed_name(name_lower: str) -> bool:
    return bool(re.search(r'_processed(\\b|\\s*\\(|_|-)', name_lower))

def _collect_image_paths(input_folder, output_folder, ignore_processed_images=True):
    in_dir  = Path(input_folder)
    out_dir = Path(output_folder).resolve()
    exts = {'.jpg', '.jpeg', '.png'}
    raw = []
    for p in in_dir.iterdir():
        if p.is_file() and p.suffix.lower() in exts:
            raw.append(p)
    seen = set()
    paths = []
    for p in raw:
        try:
            rp = p.resolve()
        except Exception:
            rp = p.absolute()
        key = str(rp).lower()
        try:
            if out_dir in rp.parents or rp == out_dir:
                continue
        except Exception:
            pass
        if ignore_processed_images and _is_processed_name(rp.stem.lower()):
            continue
        if key in seen:
            continue
        seen.add(key)
        paths.append(str(rp))
    paths.sort()
    return paths

def _perimeter_cm_from_cnt(cnt, sx, sy):
    c = cnt.reshape(-1, 2).astype(np.float64)
    c_scaled = np.empty_like(c)
    c_scaled[:,0] = c[:,0] * sx
    c_scaled[:,1] = c[:,1] * sy
    d = np.diff(np.vstack([c_scaled, c_scaled[0]]), axis=0)
    seg_lens = np.sqrt((d[:,0]**2) + (d[:,1]**2))
    return float(seg_lens.sum())

def process_images(input_folder, output_folder,
                   image_real_cm_W, image_real_cm_H,
                   lower_hsv, upper_hsv, extra_hsv,
                   margin,
                   k_open, k_close, open_iter, close_iter,
                   fill_holes,
                   min_component_area_px,
                   object_min_area_cm2,
                   rel_min_frac_of_largest,
                   max_keep,
                   outline_color_bgr,
                   outline_thickness,
                   ignore_processed_images=True):

    os.makedirs(output_folder, exist_ok=True)
    image_paths = _collect_image_paths(input_folder, output_folder,
                                      ignore_processed_images=bool(ignore_processed_images))
    rows = []

    for path in image_paths:
        filename = os.path.basename(path)
        image = cv2.imread(path)
        if image is None:
            continue

        h, w = image.shape[:2]
        area_per_pixel_cm2 = (float(image_real_cm_W) / float(w)) * (float(image_real_cm_H) / float(h))
        sx = float(image_real_cm_W) / float(w)
        sy = float(image_real_cm_H) / float(h)

        hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)

        mask = _combine_hsv_masks(hsv, lower_hsv, upper_hsv, extra_hsv)
        mask = _apply_margin(mask, margin)
        mask = _morph_open_close(mask, k_open, k_close, open_iter, close_iter)
        if bool(fill_holes):
            mask = _fill_holes(mask)

        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        # pixel area filter
        contours = [c for c in contours if cv2.contourArea(c) >= float(min_component_area_px)]
        if not contours:
            out_img = os.path.join(output_folder, os.path.splitext(filename)[0] + '_processed.jpg')
            cv2.imwrite(out_img, image)
            continue

        # candidate list with cm2
        candidates = []
        for c in contours:
            apx = float(cv2.contourArea(c))
            acm2 = apx * area_per_pixel_cm2
            candidates.append((c, apx, acm2))

        # absolute cm2 filter
        kept = [t for t in candidates if t[2] >= float(object_min_area_cm2)]

        # relative filter (optional)
        if float(rel_min_frac_of_largest) > 0 and candidates:
            largest_cm2 = max(t[2] for t in candidates)
            kept = [t for t in kept if t[2] >= float(rel_min_frac_of_largest) * largest_cm2]

        kept.sort(key=lambda t: t[2], reverse=True)
        kept = kept[:int(max_keep)]
        contours = [t[0] for t in kept]

        annotated = image.copy()
        oc = tuple(int(v) for v in outline_color_bgr)
        ot = int(outline_thickness)

        for idx, cnt in enumerate(contours, start=1):
            obj_mask = np.zeros(mask.shape, dtype=np.uint8)
            cv2.drawContours(obj_mask, [cnt], -1, 255, -1)

            area_px = float(cv2.countNonZero(obj_mask))
            area_cm2 = area_px * area_per_pixel_cm2
            perim_cm = _perimeter_cm_from_cnt(cnt, sx, sy)

            # rotated rect in CM for length/width
            pts = cnt.reshape(-1, 2).astype(np.float32)
            pts_cm = np.column_stack([pts[:,0] * sx, pts[:,1] * sy]).astype(np.float32)
            rect_cm = cv2.minAreaRect(pts_cm)
            (w_cm, h_cm) = rect_cm[1]
            length_cm = float(max(w_cm, h_cm))
            width_cm  = float(min(w_cm, h_cm))

            cv2.drawContours(annotated, [cnt], -1, oc, ot)

            rows.append({
                'File': filename,
                'Index_in_file': idx,
                'Area_cm2': round(area_cm2, 3),
                'Perimeter_cm': round(perim_cm, 3),
                'Width_cm': round(width_cm, 3),
                'Length_cm': round(length_cm, 3)
            })

        out_img = os.path.join(output_folder, os.path.splitext(filename)[0] + '_processed.jpg')
        cv2.imwrite(out_img, annotated)

    df = pd.DataFrame(rows)

    csv_path = os.path.join(output_folder, 'image_processed.csv')
    cols = ['File','Index_in_file','Area_cm2','Perimeter_cm','Width_cm','Length_cm']
    if len(df) > 0:
        df = df.reindex(columns=cols)
        df.to_csv(csv_path, index=False, encoding='utf-8-sig')
    else:
        pd.DataFrame(columns=cols).to_csv(csv_path, index=False, encoding='utf-8-sig')

    return df
"
  reticulate::py_run_string(py_code)

  ws= if (.Platform$OS.type == "windows") "\\" else "/"
  in_path= normalizePath(input_folder,  winslash = ws, mustWork = FALSE)
  out_path= normalizePath(output_folder, winslash = ws, mustWork = FALSE)

  df= reticulate::py$process_images(
    in_path,
    out_path,
    as.numeric(image_real_cm[1]), as.numeric(image_real_cm[2]),
    as.integer(lower_hsv), as.integer(upper_hsv), extra_hsv,
    as.integer(margin),
    as.integer(k_open), as.integer(k_close), as.integer(open_iter), as.integer(close_iter),
    isTRUE(fill_holes),
    as.integer(min_component_area_px),
    as.numeric(object_min_area_cm2),
    as.numeric(rel_min_frac_of_largest),
    as.integer(max_keep),
    as.integer(outline_color_bgr),
    as.integer(outline_thickness),
    isTRUE(ignore_processed_images)
  )

  reticulate::py_to_r(df)
}
# All Rights Reserved © J.K Kim (kimjk@agronomy4future.com). Last updated on 02/27/2026
