# Vision

![Vision](https://github.com/DiaaZiada/Vision/blob/master/images/Vision.jpg)

Vision is a C package for images processing stuff callable in Python. It's an assignment of [The Ancient Secrets of Computer Vision CSE 455](https://pjreddie.com/courses/computer-vision/)


**Table Content**

 * [Basic Level Of Image Processing](#basic-level-of-image-processing)
	 * [Getting and setting pixels](#getting-and-setting-pixels)
	 * [Copy image](#copy-image)
	 * [ Grayscale image](#grayscale-image)
	 * [ Shifting the image colors](#shifting-the-image-colors)
	 * [Clamping the image values](#clamping-the-image-values)
	 * [RGB and HSV](#rgb-and-hsv)
 * [Middle Level Of Image Processing](#middle-level-of-image-processing)
	 * [Image resizing](#image-resizing) 
	 * [Image filtering with convolutions](#image-filtering-with-convolutions) 
		* [Blurring](#blurring) 
		* [Blur and Resizing](#blur-and-resizing)
		* [Gaussian](#gaussian)
		* [Frequency and Reconstruction](#frequency-and-reconstruction)
	* [Image Features](#image-features)
		* [Sobel filters](#sobel-filters) 
		* [Colorized Representation](#colorized-representation)
 * [Panorama](#panorama)
	 * [Harris detector](#harris-detector) 
	 * [Patch matching](#patch-matching)
	 * [Combine the images with a homography](#combine-the-images-with-a-homography)
 4. Optical Flow
 5. Neural Network
	
Original Image

![Original Image](https://github.com/DiaaZiada/Vision/blob/master/data/dog.jpg)

## Basic Level Of Image Processing

### Getting and setting pixels
**Getting pixel value**
```python
from vision import *
im = load_image("data/dog.jpg")
pixel = get_pixel(im, 0, 0, 0)
print(pixel)
```
```python
=> 0.0470588244497776    
```
**Setting pixel to a value**

```python
from vision import *
# remove all red colors from the image
im = load_image("data/dog.jpg")
for row in range(im.h):
    for col in range(im.w):
        set_pixel(im, col, row, 0, 0)
save_image(im, "output/no_red_dog")
```

![no red](https://github.com/DiaaZiada/Vision/blob/master/output/no_red_dog.jpg)

### Copy image
```python
from vision import *
im = load_image("data/dog.jpg")
im2 = copy_image(im) 
```
### Grayscale image
```python
from vision import *
im = load_image("data/dog.jpg")
graybar = rgb_to_grayscale(im)
save_image(graybar, "output/gray_dog")
```
![gray](https://github.com/DiaaZiada/Vision/blob/master/output/gray_dog.jpg)

### Shifting the image colors
```python
from vision import *
im = load_image("data/dog.jpg")
shift_image(im, 0, .4)
shift_image(im, 1, .4)
shift_image(im, 2, .4)
save_image(im, "output/shifted_dog")
```
![shifted](https://github.com/DiaaZiada/Vision/blob/master/output/shifted_dog.jpg)

### Clamping the image values
```python
from vision import *
im = load_image("data/dog.jpg")
shift_image(im, 0, .4)
shift_image(im, 1, .4)
shift_image(im, 2, .4)
clamp_image(im)
save_image(im, "output/ligth_fixed_dog")
```
![Clamping](https://github.com/DiaaZiada/Vision/blob/master/output/light_fixed_dog.jpg)

### RGB and HSV
```python
from vision import *
im = load_image("data/dog.jpg")
rgb_to_hsv(im)
shift_image(im, 1, .2)
clamp_image(im)
hsv_to_rgb(im)
save_image(im, "output/rgb_hsv_rgb_dog")
```
![shifted](https://github.com/DiaaZiada/Vision/blob/master/output/rgb_hsv_rgb_dog.jpg)

## Middle Level Of Image Processing

### Image resizing
#### Maximizing

**Nearest Neighbor** 

```python
from vision import *
im = load_image("data/dog.jpg")
a = nn_resize(im, im.w*4, im.h*4)
save_image(a, "output/4x_nn_dog")
```
![maximizing nn ](https://github.com/DiaaZiada/Vision/blob/master/output/4x_nn_dog.jpg)

**Bilinear**
```python
from vision import *
im = load_image("data/dog.jpg")
a = bilinear_resize(im, im.w*4, im.h*4)
save_image(a, "output/4x_bl_dog")
```
![maximizing bl ](https://github.com/DiaaZiada/Vision/blob/master/output/4x-bl_dog.jpg)


#### Minimizing 

**Nearest Neighbor** 

```python
from vision import *
im = load_image("data/dog.jpg")
a = nn_resize(im, im.w//7, im.h//7)
save_image(a, "output/7th-nn_dog")
```
![minimizing nn ](https://github.com/DiaaZiada/Vision/blob/master/output/7th-nn_dog.jpg)

**Bilinear**
```python
from vision import *
im = load_image("data/dog.jpg")
a = bilinear_resize(im, im.w//7, im.h//7)
save_image(a, "output/7th-bl_dog")
```
![minimizing bl ](https://github.com/DiaaZiada/Vision/blob/master/output/7th-bl_dog.jpg)

### Image filtering with convolutions

#### Blurring 
```python
from vision import *
im = load_image("data/dog.jpg")
f = make_box_filter(7)
blur = convolve_image(im, f, 1)
save_image(blur, "output/box7_dog")
```
![Blurring](https://github.com/DiaaZiada/Vision/blob/master/output/box7_dog.jpg)

#### Blur and Resizing 

```python
from vision import *
im = load_image("data/dog.jpg")
f = make_box_filter(7)
blur = convolve_image(im, f, 1)
thumb = nn_resize(blur, blur.w//7, blur.h//7)
save_image(thumb, "output/thumb_dog")
```
![ Blur and Resizing ](https://github.com/DiaaZiada/Vision/blob/master/output/thumb_dog.jpg)

#### Gaussian
```python
from vision import *
im = load_image("data/dog.jpg")
f = make_gaussian_filter(2)
blur = convolve_image(im, f, 1)
save_image(blur, "output/gauss2_dog")
```
![Gaussian](https://github.com/DiaaZiada/Vision/blob/master/output/gauss2_dog.jpg)

#### Frequency and Reconstruction
```python
from vision import *
im = load_image("data/dog.jpg")
f = make_gaussian_filter(2)
lfreq = convolve_image(im, f, 1)
hfreq = im - lfreq
reconstruct = lfreq + hfreq
save_image(lfreq, "output/low_frequency_dog.jpg")
save_image(hfreq, "output/high_frequency_dog.jpg")
save_image(reconstruct, "output/reconstruct_dog.jpg")
```
low frequency

![low frequency](https://github.com/DiaaZiada/Vision/blob/master/output/low_frequency_dog.jpg)

high frequency

![high frequency ](https://github.com/DiaaZiada/Vision/blob/master/output/high_frequency_dog.jpg)

reconstruct dog

![reconstruct dog](https://github.com/DiaaZiada/Vision/blob/master/output/reconstruct_dog.jpg)

### Image Features
#### Sobel filters
```python
from vision import *
im = load_image("data/dog.jpg")
res = sobel_image(im)
mag = res[0]
feature_normalize(mag)
save_image(mag, "output/magnitude")
```
![magnitude](https://github.com/DiaaZiada/Vision/blob/master/output/magnitude.jpg)

#### Colorized Representation
```python
from vision import *
im = load_image("data/dog.jpg")
res = sobel_image(im)
mag = res[0]
feature_normalize(mag)
img = colorize_sobel(im)
save_image(img, "output/colored_magnitude_dog")
```
![colored magnitude](https://github.com/DiaaZiada/Vision/blob/master/output/colored_magnitude_dog.jpg)

## Panorama
### Harris detector
```python
from vision import *
im = load_image("data/Rainier1.png")
detect_and_draw_corners(im, 2, 50, 3)
save_image(im, "output/corners")
```
![harris detector](https://github.com/DiaaZiada/Vision/blob/master/output/corners.jpg)

### Patch matching
```python
from vision import *
a = load_image("data/Rainier1.png")
b = load_image("data/Rainier2.png")
m = find_and_draw_matches(a, b, 2, 50, 3)
save_image(m, "output/matches")
```
![patch matching](https://github.com/DiaaZiada/Vision/blob/master/output/matches.jpg)

### Combine the images with a homography
```python
from vision import *
im1 = load_image("data/Rainier1.png")
im2 = load_image("data/Rainier2.png")
pan = panorama_image(im1, im2, thresh=50)
save_image(pan, "output/easy_panorama")
```
![combine the images with a homography](https://github.com/DiaaZiada/Vision/blob/master/output/easy_panorama.jpg)

**Final panorama image**
![final panorama image](https://github.com/DiaaZiada/Vision/blob/master/output/rainier_panoram.jpg)







