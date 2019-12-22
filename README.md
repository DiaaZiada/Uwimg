# Vision

![Vision](https://github.com/DiaaZiada/Vision/blob/master/images/Vision.jpg)

Vision is C package for images processing stuff callable in Python

**Table Content**

 

 * [Basic Level Of Image Processing](#basic-level-of-image-processing)
	 * [Getting and setting pixels](#getting-and-setting-pixels)
 2. Middle Level Of Image Processing]
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
# remove all red colors from the image
im = load_image("data/dog.jpg")
for row in range(im.h):
    for col in range(im.w):
        set_pixel(im, col, row, 0, 0)
save_image(im, "output/dog")
=>
```

![no red](https://github.com/DiaaZiada/Vision/blob/master/output/dog.jpg)




 
