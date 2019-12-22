#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO Fill in
    return get_pixel(im, round(x), round(y), c);   
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix that first line)
    image new_image = make_image(w, h, im.c);
    float old_top_x, old_top_y, old_bottom_x, old_bottom_y ,
    new_top_x, new_top_y, new_bottom_x, new_bottom_y, ax, bx,
    ay, by, new_x, new_y;
    
    old_top_x = -0.5;
    old_top_y = -0.5;
    old_bottom_x = im.w - 0.5;
    old_bottom_y = im.h - 0.5;

    new_top_x = -0.5;
    new_top_y = -0.5;
    new_bottom_x = w - 0.5;
    new_bottom_y = h - 0.5;


    ax = (old_bottom_x - old_top_x) / (new_bottom_x - new_top_x);
    ay = (old_bottom_y - old_top_y) / (new_bottom_y - new_top_y);
    
    bx = old_top_x - (ax * new_top_x);
    by = old_top_y - (ay * new_top_y);

    for(int c=0; c < new_image.c; c++)
        for(int y=0; y<new_image.h; y++)
            for(int x=0; x<new_image.w; x++){
                new_x = ax*x + bx;
                new_y = ay*y + by;
                set_pixel(new_image, x, y, c, nn_interpolate(im, new_x, new_y, c));
            }


    return new_image;


}

float bilinear_interpolate(image im, float x, float y, int c)
{
    int right = (int) ceilf(x);
    int left = (int) floorf(x);
    int top = (int) floorf(y);
    int buttom = (int) ceilf(y);

    float q1 = (buttom - y) * get_pixel(im,left,top,c) + (y - top) * get_pixel(im,left,buttom,c);

    float q2 = (buttom - y) * get_pixel(im,right,top,c) + (y - top) * get_pixel(im,right,buttom,c);

    return (right - x) * q1 + (x - left) * q2;
}


image bilinear_resize(image im, int w, int h)
{
    // TODO
    image new_image = make_image(w, h, im.c);
    float old_top_x, old_top_y, old_bottom_x, old_bottom_y ,
    new_top_x, new_top_y, new_bottom_x, new_bottom_y, ax, bx,
    ay, by, new_x, new_y;
    
    old_top_x = -0.5;
    old_top_y = -0.5;
    old_bottom_x = im.w - 0.5;
    old_bottom_y = im.h - 0.5;

    new_top_x = -0.5;
    new_top_y = -0.5;
    new_bottom_x = w - 0.5;
    new_bottom_y = h - 0.5;


    ax = (old_bottom_x - old_top_x) / (new_bottom_x - new_top_x);
    ay = (old_bottom_y - old_top_y) / (new_bottom_y - new_top_y);
    
    bx = old_top_x - (ax * new_top_x);
    by = old_top_y - (ay * new_top_y);

    for(int c=0; c < new_image.c; c++)
        for(int y=0; y<new_image.h; y++)
            for(int x=0; x<new_image.w; x++){
                new_x = ax*x + bx;
                new_y = ay*y + by;
                set_pixel(new_image, x, y, c, bilinear_interpolate(im, new_x, new_y, c));
            }


    return new_image;


}

