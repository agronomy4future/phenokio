# ============================================================
# Package: phenokio
# Title: Automated Grain Morphometric Analysis and Seed Phenotyping
# Description: High-throughput grain analysis using R and OpenCV.
# All Rights Reserved © J.K Kim (kimjk@agronomy4future.com)
# Last updated on 02/27/2026
# ============================================================

if (!requireNamespace("reticulate", quietly = TRUE)) install.packages("reticulate")
library(reticulate)

#' Internal helper to ensure Python dependencies are installed
ensure_pydeps <- function() {
  invisible(reticulate::py_available(initialize = TRUE))
  for (pkg in c("numpy", "pandas")) {
    if (!reticulate::py_module_available(pkg)) reticulate::py_install(pkg, pip = TRUE)
  }
  if (!reticulate::py_module_available("cv2")) reticulate::py_install("opencv-python-headless", pip = TRUE)
}

#' Resolve any R color name or BGR vector to OpenCV-compatible BGR
.resolve_color_bgr <- function(col) {
  # 1. If numeric vector c(B, G, R) is provided
  if (is.numeric(col) && length(col) == 3) return(as.integer(col))

  # 2. Use R's internal color database to support all names (orange, pink, etc.)
  rgb_val <- tryCatch(grDevices::col2rgb(col), error = function(e) NULL)

  if (is.null(rgb_val)) {
    message("Unknown color name. Defaulting to green.")
    return(c(0L, 255L, 0L))
  }

  # Convert RGB to BGR for OpenCV
  return(as.integer(rev(as.vector(rgb_val))))
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
phenokio <- function(
    input_folder,
    output_folder,
    image_real_cm = c(20, 20),

    # --- segmentation (HSV) ---
    lower_hsv = c(0L, 0L, 0L),
    upper_hsv = c(180L, 255L, 255L),
    extra_hsv = NULL,

    # --- border exclusion ---
    margin = 0L,

    # --- morphology ---
    k_open  = c(5L, 5L),
    k_close = c(7L, 7L),
    open_iter  = 1L,
    close_iter = 2L,
    fill_holes = FALSE,

    # --- filtering (Noise & Line prevention) ---
    min_component_area_px = 200L,
    object_min_area_cm2 = 0.01,
    rel_min_frac_of_largest = 0.00,
    max_aspect_ratio = 4.0,        # [New] If length/width > 4.0, it's likely a line, not a seed
    max_keep = 2000L,

    # --- output overlay style ---
    outline_color = "green",       # Now supports "orange", "purple", "skyblue", etc.
    outline_thickness = 2L,

    ignore_processed_images = TRUE
) {
  ensure_pydeps()

  if (length(image_real_cm) == 1) image_real_cm <- rep(image_real_cm, 2)
  dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

  # Support for ALL R color names
  outline_color_bgr <- .resolve_color_bgr(outline_color)

  py_code <- "
import cv2
import numpy as np
import os, pandas as pd
from pathlib import Path
import re

def _as_uint8_triplet(x):
    arr = np.array(x, dtype=np.int64).reshape(-1)
    return np.clip(arr, 0, 255).astype(np.uint8)

def _combine_hsv_masks(hsv_img, lower_hsv, upper_hsv, extra_hsv):
    base = cv2.inRange(hsv_img, _as_uint8_triplet(lower_hsv), _as_uint8_triplet(upper_hsv))
    if extra_hsv is None: return base
    masks = [base]
    try: iterable = list(extra_hsv)
    except: iterable = []
    for item in iterable:
        try:
            lo = item.get('lower', item[0]); hi = item.get('upper', item[1])
            masks.append(cv2.inRange(hsv_img, _as_uint8_triplet(lo), _as_uint8_triplet(hi)))
        except: continue
    out = masks[0]
    for m in masks[1:]: out = cv2.bitwise_or(out, m)
    return out

def _apply_margin(mask, margin):
    m = int(margin)
    if m <= 0: return mask
    h, w = mask.shape[:2]
    mask[:m, :] = 0; mask[h-m:h, :] = 0
    mask[:, :m] = 0; mask[:, w-m:w] = 0
    return mask

def _morph_open_close(mask, k_open, k_close, open_iter, close_iter):
    if k_open[0] > 0:
        mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, np.ones(tuple(k_open), np.uint8), iterations=int(open_iter))
    if k_close[0] > 0:
        mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, np.ones(tuple(k_close), np.uint8), iterations=int(close_iter))
    return mask

def _fill_holes(mask):
    h, w = mask.shape[:2]
    flood = cv2.bitwise_not(mask)
    ffmask = np.zeros((h+2, w+2), dtype=np.uint8)
    cv2.floodFill(flood, ffmask, (0,0), 255)
    return cv2.bitwise_or(mask, cv2.bitwise_not(flood))

def process_images(input_folder, output_folder, image_real_cm_W, image_real_cm_H,
                   lower_hsv, upper_hsv, extra_hsv, margin, k_open, k_close,
                   open_iter, close_iter, fill_holes, min_px, min_cm2, rel_frac,
                   max_ar, max_keep, oc_bgr, ot, ignore_proc):

    exts = {'.jpg', '.jpeg', '.png'}
    paths = [str(p) for p in Path(input_folder).iterdir() if p.suffix.lower() in exts]
    rows = []

    for path in paths:
        filename = os.path.basename(path)
        if ignore_proc and '_processed' in filename.lower(): continue

        img = cv2.imread(path)
        if img is None: continue
        h, w = img.shape[:2]
        sx, sy = float(image_real_cm_W)/w, float(image_real_cm_H)/h
        area_factor = sx * sy

        hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
        mask = _combine_hsv_masks(hsv, lower_hsv, upper_hsv, extra_hsv)
        mask = _apply_margin(mask, margin)
        mask = _morph_open_close(mask, k_open, k_close, open_iter, close_iter)
        if fill_holes: mask = _fill_holes(mask)

        cnts, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        candidates = []
        for c in cnts:
            px_area = cv2.contourArea(c)
            cm2_area = px_area * area_factor
            if px_area < min_px or cm2_area < min_cm2: continue

            # Line filter: Aspect Ratio (Length/Width)
            pts_cm = (c.reshape(-1, 2) * [sx, sy]).astype(np.float32)
            rect = cv2.minAreaRect(pts_cm)
            dim1, dim2 = rect[1]
            if dim1 == 0 or dim2 == 0: continue
            ar = max(dim1, dim2) / min(dim1, dim2)

            if ar > max_ar: continue # Skip if it looks like a line
            candidates.append((c, cm2_area, max(dim1, dim2), min(dim1, dim2)))

        if rel_frac > 0 and candidates:
            max_area = max(t[1] for t in candidates)
            candidates = [t for t in candidates if t[1] >= rel_frac * max_area]

        candidates.sort(key=lambda x: x[1], reverse=True)
        candidates = candidates[:max_keep]

        annotated = img.copy()
        for idx, (cnt, area, l_cm, w_cm) in enumerate(candidates, 1):
            cv2.drawContours(annotated, [cnt], -1, tuple(int(v) for v in oc_bgr), int(ot))
            rows.append({'File': filename, 'Index': idx, 'Area_cm2': round(area, 3),
                         'Width_cm': round(w_cm, 3), 'Length_cm': round(l_cm, 3)})

        cv2.imwrite(os.path.join(output_folder, os.path.splitext(filename)[0] + '_processed.jpg'), annotated)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(output_folder, 'image_processed.csv'), index=False, encoding='utf-8-sig')
    return df
"
  reticulate::py_run_string(py_code)

  df <- reticulate::py$process_images(
    normalizePath(input_folder, winslash = "/", mustWork = FALSE),
    normalizePath(output_folder, winslash = "/", mustWork = FALSE),
    as.numeric(image_real_cm[1]), as.numeric(image_real_cm[2]),
    as.integer(lower_hsv), as.integer(upper_hsv), extra_hsv,
    as.integer(margin), as.integer(k_open), as.integer(k_close),
    as.integer(open_iter), as.integer(close_iter), isTRUE(fill_holes),
    as.integer(min_component_area_px), as.numeric(object_min_area_cm2),
    as.numeric(rel_min_frac_of_largest), as.numeric(max_aspect_ratio),
    as.integer(max_keep), as.integer(outline_color_bgr),
    as.integer(outline_thickness), isTRUE(ignore_processed_images)
  )

  return(reticulate::py_to_r(df))
}
# All Rights Reserved © J.K Kim (kimjk@agronomy4future.com). Last updated on 02/27/2026
