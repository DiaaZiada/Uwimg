# Vision

![Vision](https://github.com/DiaaZiada/Vision/blob/master/images/Vision.jpg)

Vision is C package for images processing stuff callable in Python

**Table Content**

 * [Basic Level Of Image Processing](#basic-level-of-image-processing)
	 * [Getting and setting pixels](#getting-and-setting-pixels)
	 * [Copy image](#copy-image)
	 * [ Grayscale image](#grayscale-image)
	 * [ Shifting the image colors](#shifting-the-image-colors)
	 * [Clamping the image values](#clamping-the-image-values)
	 * [RGB and HSV](#rgb-and-hsv)
 * [Middle Level Of Image Processing](#middle-level-of-image-processing)
 3. Panorama
 4. Optical Flow
 5. Neural Network
	
## Basic Level Of Image Processing

### Getting and setting pixels
**Getting pixel value**
```python
from vision import *
im = load_image("data/dog.jpg")
pixel = get_pixel(im, 0, 0, 0)
print(pixel)
=> 0.0470588244497776    
```
**Setting pixel to a value**

remove all red colors from the image
```python
from vision import *
# remove all red colors from the image
im = load_image("data/dog.jpg")
for row in range(im.h):
    for col in range(im.w):
        set_pixel(im, col, row, 0, 0)
save_image(im, "output/no_red_dog")
=>
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
=>
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
=>
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
=>
```
![shifted](https://github.com/DiaaZiada/Vision/blob/master/output/ligth_fixed_dog.jpg)

### RGB and HSV
```python
from vision import *
im = load_image("data/dog.jpg")
rgb_to_hsv(im)
shift_image(im, 1, .2)
clamp_image(im)
hsv_to_rgb(im)
save_image(im, "output/rgb_hsv_rgb_dog")
=>
```
![shifted](https://github.com/DiaaZiada/Vision/blob/master/output/rgb_hsv_rgb_dog.jpg)

## Middle Level Of Image Processing

###
 
