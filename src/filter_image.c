#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define M_PI 3.14159265358979323846

void l1_normalize(image im)
{
    /**
    Function: apply lq normalization to an image
    Args:
        im: struct contains images shape(weight, height, and channels)(im.w, im.h, im.c) and image pixels (im.data)  
    */
    for(int x=0; x<im.w; x++)
        for(int y=0; y<im.h; y++)
            for(int c=0; c<im.c; c++)
                set_pixel(im,x,y,c,get_pixel(im,x,y,c)/(1.0*im.w*im.h));
}

image make_box_filter(int w)
{   
    /**
    Function: retruns an filter with dims (WxW)
    Args:
        w: integer for the dim of the filter
    Return:
        filter: struct contains filer shape(weight, height, and 1)(filer.w, filer.h, filer.c) and filter  pixels (filer.data)  
    */
    image filter = make_image(w,w,1);
    for(int x=0; x<filter.w; x++)
        for(int y=0; y<filter.h; y++)
            for(int c=0; c<filter.c; c++)
                set_pixel(filter,x,y,c,1.0);
    l1_normalize(filter);
    return filter;
}
image convolve_image(image im, image filter, int preserve)
{   
    /**
    Function: apply convolution operation
    Args:
        im: struct contains images shape(weight, height, and channels)(im.w, im.h, im.c) and image pixels (im.data)
        filter: struct contains filer shape(weight, height, and channels)(filer.w, filer.h, filer.c) and filter  pixels (filer.data)  
        preserve: integer if 1 then keep input chennels to the output else output channels be 1
     Return:
        new_image: (output of the operation) struct contains images shape(weight, height, and channels)(new_image.w, new_image.h, new_image.c) and image pixels (new_image.data)
    
    */
    assert(filter.c == im.c || filter.c == 1);
    image new_image = make_image(im.w, im.h, im.c);
    for (int x = 0; x < im.w; x++) {
        for (int y = 0; y < im.h; y++) {
            for (int c = 0; c < im.c; c++) {
                float q = 0;
                for (int filter_x = 0; filter_x < filter.w; filter_x++) {
                    for (int filter_y = 0; filter_y < filter.h; filter_y++) {
                        float value = filter.data[filter_x + filter_y*filter.w + filter.w*filter.h*(filter.c-1)*c];
                        int fx = x - filter.w / 2 + filter_x;
                        int fy = y - filter.h / 2 + filter_y;
                        q += get_pixel(im, fx, fy, c) * value;
                    }
                }
                new_image.data[x + y*im.w + c*im.w*im.h] = q;
            }
        }
    }
    if (!preserve) {
        image new_image_ = make_image(im.w, im.h, 1);
        for (int w = 0; w < im.w; w++) {
            for (int h = 0; h < im.h; h++) {
                float q = 0;
                for (int c = 0; c < im.c; c++) {
                q += new_image.data[w + h*im.w + im.w*im.h*c];
                }
                new_image_.data[w + h*im.w] = q;
            }
        }
        return new_image_;
    }
    else {
        return new_image;
    }
}


image make_highpass_filter()
{
    /**
    Function: retruns highpass filter
    Return:
        filter: struct contains filer shape(weight, height, and 1)(filer.w, filer.h, filer.c) and filter  pixels (filer.data)  
    */
    image f = make_box_filter(3);
    set_pixel(f,0,0,0,0);
    set_pixel(f,1,0,0,-1);
    set_pixel(f,2,0,0,0);
    set_pixel(f,0,1,0,-1);
    set_pixel(f,1,1,0,4);
    set_pixel(f,2,1,0,-1);
    set_pixel(f,0,2,0,0);
    set_pixel(f,1,2,0,-1);
    set_pixel(f,2,2,0,0);
    return f;
}

image make_sharpen_filter()
{   /**
    Function: retruns sharpen filter
    Return:
        filter: struct contains filer shape(weight, height, and 1)(filer.w, filer.h, filer.c) and filter  pixels (filer.data)  
    */
    image f = make_box_filter(3);
    set_pixel(f,0,0,0,0);
    set_pixel(f,1,0,0,-1);
    set_pixel(f,2,0,0,0);
    set_pixel(f,0,1,0,-1);
    set_pixel(f,1,1,0,5);
    set_pixel(f,2,1,0,-1);
    set_pixel(f,0,2,0,0);
    set_pixel(f,1,2,0,-1);
    set_pixel(f,2,2,0,0);
    return f;
}

image make_emboss_filter()
{   /**
    Function: retruns emboss filter
    Return:
        filter: struct contains filer shape(weight, height, and 1)(filer.w, filer.h, filer.c) and filter  pixels (filer.data)  
    */
    image f = make_box_filter(3);
    set_pixel(f,0,0,0,-2);
    set_pixel(f,1,0,0,-1);
    set_pixel(f,2,0,0,0);
    set_pixel(f,0,1,0,-1);
    set_pixel(f,1,1,0,1);
    set_pixel(f,2,1,0,1);
    set_pixel(f,0,2,0,0);
    set_pixel(f,1,2,0,1);
    set_pixel(f,2,2,0,2);
    return f;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    /**
    Function: retruns gaussian filter
    Args:
        sigma: float number of sigma value
    Return:
        f: struct contains filer shape(weight, height, and 1)(f.w, f.h, f.c) and filter  pixels (f.data)  
    */
    image f = make_box_filter(sigma*6+1);
    float g,gx,gy;
    for (int x=0; x<f.w; x++)
        for (int y=0; y<f.h; y++){
            gx = x - 3 * sigma;
            gy = y - 3 * sigma;
            g = (1./(2*M_PI* powf(sigma,2))) * powf(M_E,-(powf(gx,2)+powf(gy,2))/(2*powf(sigma,2)));
            set_pixel(f,x,y,0,g);
   
        }

    return f;
}

image add_image(image a, image b)
{   /**
    Function: add two images' pixel 
    Args:
        a: (image one) struct contains images shape(weight, height, and channels)(a.w, a.h, a.c) and image pixels (a.data)  
        b: (image two)struct contains images shape(weight, height, and channels)(b.w, b.h, b.c) and image pixels (b.data)  
    Return:
        new_image: (produced image)struct contains images shape(weight, height, and channels)(new_image.w, new_image.h, new_image.c) and image pixels (new_image.data) 
    */
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    image new_image = make_image(a.w,a.h,a.c);
    float new_pix=0.;
    for (int x = 0; x < new_image.w; x++) 
        for (int y = 0; y < new_image.h; y++) 
            for (int c = 0; c < new_image.c; c++) {
                new_pix = get_pixel(a,x,y,c) + get_pixel(b,x,y,c);
                set_pixel(new_image, x, y, c, new_pix);
            }

    return new_image;
}

image sub_image(image a, image b)
{
    /**
    Function: subtract two images' pixel
    Args:
        a: (image one) struct contains images shape(weight, height, and channels)(a.w, a.h, a.c) and image pixels (a.data)  
        b: (image two)struct contains images shape(weight, height, and channels)(b.w, b.h, b.c) and image pixels (b.data)  
    Return:
        new_image: (produced image)struct contains images shape(weight, height, and channels)(new_image.w, new_image.h, new_image.c) and image pixels (new_image.data) 
    */

    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    image new_image = make_image(a.w,a.h,a.c);
    float new_pix=0.;
    for (int x = 0; x < new_image.w; x++) 
        for (int y = 0; y < new_image.h; y++) 
            for (int c = 0; c < new_image.c; c++) {
                new_pix = get_pixel(a,x,y,c) - get_pixel(b,x,y,c);
                set_pixel(new_image, x, y, c, new_pix);
            }

    return new_image;
}

image make_gx_filter()
{   /**
    Function: retruns sobel x filter
    Return:
        filter: struct contains filer shape(weight, height, and 1)(filer.w, filer.h, filer.c) and filter  pixels (filer.data)  
    */
    image f = make_box_filter(3);
    set_pixel(f,0,0,0,-1);
    set_pixel(f,1,0,0,0);
    set_pixel(f,2,0,0,1);
    set_pixel(f,0,1,0,-2);
    set_pixel(f,1,1,0,0);
    set_pixel(f,2,1,0,2);
    set_pixel(f,0,2,0,-1);
    set_pixel(f,1,2,0,0);
    set_pixel(f,2,2,0,1);
    return f;
}

image make_gy_filter()
{   /**
    Function: retruns sobel y filter
    Return:
        filter: struct contains filer shape(weight, height, and 1)(filer.w, filer.h, filer.c) and filter  pixels (filer.data)  
    */
    image f = make_box_filter(3);
    set_pixel(f,0,0,0,-1);
    set_pixel(f,1,0,0,-2);
    set_pixel(f,2,0,0,-1);
    set_pixel(f,0,1,0,0);
    set_pixel(f,1,1,0,0);
    set_pixel(f,2,1,0,0);
    set_pixel(f,0,2,0,1);
    set_pixel(f,1,2,0,2);
    set_pixel(f,2,2,0,1);
    return f;
}

void feature_normalize(image im)
{   /**
    Function: apply normalization to an image
    Args:
        im: struct contains images shape(weight, height, and channels)(im.w, im.h, im.c) and image pixels (im.data)
    */
    float max, min;
    max= -INFINITY;
    min= INFINITY;
    for (int i=0; i<im.w*im.h*im.c; i++){
        if(im.data[i]>max)
            max = im.data[i];
        if(im.data[i]<min)
            min = im.data[i];
    }
    if(max - min !=0) 
        for (int i=0; i<im.w*im.h*im.c; i++)
            im.data[i] = (im.data[i] - min) / (max - min);
    else
         for (int i=0; i<im.w*im.h*im.c; i++)
            im.data[i] = 0.;

}

image *sobel_image(image im)
{
    /**
    Function: apply sobel operation to an image
    Args:
        im: struct contains images shape(weight, height, and channels)(im.w, im.h, im.c) and image pixels (im.data)
    Return:
        sobel: array of grayscale images 1st contains magnitude, 2nd contains direction
    */
    float magnitude, direction;
    image gx = make_gx_filter();
    image gy = make_gy_filter();
    
    image * sobel = calloc(2, sizeof(image));

    image sobel_x = convolve_image(im, gx , 0);
    image sobel_y = convolve_image(im, gy , 0);
    
    sobel[0] = make_image(im.w,im.h,1);
    sobel[1] = make_image(im.w,im.h,1);

    for(int x=0; x<im.w; x++)
        for(int y=0; y<im.h; y++){
            magnitude = sqrtf(powf(get_pixel(sobel_x,x,y,0),2) + powf(get_pixel(sobel_y,x,y,0),2));
            direction =  atan2f(get_pixel(sobel_y,x,y,0), get_pixel(sobel_x,x,y,0));
            set_pixel(sobel[0],x,y,0,magnitude);
            set_pixel(sobel[1],x,y,0,direction);
        }   
    return sobel;
}

image colorize_sobel(image im)
{
    /**
    Function: colorize the sobel image data
    Args:
        im: struct contains images shape(weight, height, and channels)(im.w, im.h, im.c) and image pixels (im.data)
    Return:
        (colored image) struct contains images shape(weight, height, and channels)(im.w, im.h, im.c) and image pixels (im.data)
    */
    image * sobel = sobel_image(im);
    feature_normalize(sobel[0]);
    feature_normalize(sobel[1]);

    image colored_sobel = make_image(im.w, im.h, 3);
    for(int c=0; c<colored_sobel.c; c++)
        for(int y=0; y<colored_sobel.h; y++)
            for(int x=0; x<colored_sobel.w; x++){
                if(c==0)
                    set_pixel(colored_sobel, x, y, c, get_pixel(sobel[1],x,y,0));
                else
                    set_pixel(colored_sobel, x, y, c, get_pixel(sobel[0],x,y,0));
            }

    hsv_to_rgb(colored_sobel);
    
    return convolve_image(colored_sobel, make_gaussian_filter(1), 1);
}
