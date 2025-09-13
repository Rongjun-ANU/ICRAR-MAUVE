'''
APPROACH:
This python code performs color-based segmentation using HSV color space, focusing on the Hue component
with additional filtering on Saturation and Value components.

The segmentation process follows these steps:
1. Color Space Conversion: Convert RGB image to HSV color space, which separates color information (Hue) 
   from illumination intensity (Value) and color purity (Saturation).
2. Thresholding: Create a binary mask where pixels with Hue values between the lower and upper thresholds 
   (and with Saturation and Value ≥ 100) are set to white (255) and the rest are set to black (0).
3. Mask Application: Apply the binary mask to the original image to extract only the color regions of interest.
4. Binarization: Convert the segmented image to a binary image using thresholding.
5. Erosion: Apply morphological erosion to remove small details from the foreground.
6. Morphological Closing: Apply morphological closing to remove small holes in the foreground objects.
7. Ellipse Extraction: Find contours and fit an ellipse to the largest contour to extract the galaxy shape.
'''

# Import necessary libraries
import cv2
import numpy as np

def segment_image(image_path, lower_hue, upper_hue):
    img = cv2.imread(image_path)
    if img is None:
        print(f"Error: Could not read image at {image_path}")
        return None, None
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    mask = cv2.inRange(hsv,
                       np.array([lower_hue, 100, 100]),
                       np.array([upper_hue, 255, 255]))
    segmented = cv2.bitwise_and(img, img, mask=mask)
    return segmented, mask

def binarize_image(image):
    if len(image.shape) == 3:
        gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    else:
        gray = image
    _, binary = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY)
    return binary

def extract_ellipse(binary_image):
    contours, _ = cv2.findContours(binary_image,
                                   cv2.RETR_EXTERNAL,
                                   cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        print("No contours found in the image")
        return None, None
    largest = max(contours, key=cv2.contourArea)
    if len(largest) < 5:
        print("Not enough points to fit an ellipse")
        return None, None

    ellipse = cv2.fitEllipse(largest)
    ellipse_image = np.zeros_like(binary_image)
    cv2.ellipse(ellipse_image, ellipse, 255, 2)
    (cx, cy), (a, b), angle = ellipse
    print(f"Center: ({cx:.2f}, {cy:.2f})")
    print(f"Axes: {max(a,b):.2f}, {min(a,b):.2f}")
    print(f"Angle: {angle:.2f}")
    print(f"Eccentricity: {(1 - (min(a,b)/max(a,b))**2)**0.5:.4f}")
    return ellipse, ellipse_image

if __name__ == "__main__":
    image_path = "ic3392_g.png"
    lower_hue, upper_hue = 83, 108

    segmented_image, mask = segment_image(image_path, lower_hue, upper_hue)
    if segmented_image is None:
        exit(1)

    # Binarize the mask
    binary_image = binarize_image(mask)

    # Draw all contours and highlight the largest
    contours, _ = cv2.findContours(binary_image,
                                   cv2.RETR_EXTERNAL,
                                   cv2.CHAIN_APPROX_SIMPLE)
    contour_vis = np.zeros_like(segmented_image)
    cv2.drawContours(contour_vis, contours, -1, (0, 0, 255), 2)
    if contours:
        largest = max(contours, key=cv2.contourArea)
        cv2.drawContours(contour_vis, [largest], -1, (255, 0, 255), 3)

    # Fit ellipse directly on the largest contour
    ellipse, ellipse_image = extract_ellipse(binary_image)

    # Save and show results
    cv2.imwrite("segmented_" + image_path, segmented_image)
    cv2.imshow("Segmented Image", segmented_image)
    cv2.imwrite("binary_" + image_path, binary_image)
    cv2.imshow("Binary Image", binary_image)
    cv2.imwrite("contours_" + image_path, contour_vis)
    cv2.imshow("Contours", contour_vis)
    if ellipse_image is not None:
        cv2.imwrite("ellipse_" + image_path, ellipse_image)
        orig = cv2.imread(image_path)
        cv2.ellipse(orig, ellipse, (0, 255, 0), 2)
        cv2.imwrite("ellipse_fit_" + image_path, orig)
        cv2.imshow("Ellipse Image", ellipse_image)
        cv2.imshow("Ellipse Fit", orig)
        print("Ellipse image saved.")
    else:
        print("No ellipse could be fitted to the contours.")

    print("Processed images saved.")
    cv2.waitKey(0)
    cv2.destroyAllWindows()
