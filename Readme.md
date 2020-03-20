Fourier Epicircle Animation for 2d Curves
====



## Mathematical Background

A 2d  curve is $t\mapsto [x(t),y(t)],\ t\in[0,1]$ , essentially we can use complex representation $z(t)=x(t)+iy(t)$. 

A closed curve is a 2d periodic signal, thus we can do Fourier analysis on it. Analytically, we want a decomposition 
$$
z(t)=\sum_{k=0} U_ke^{-i2\pi kt}
$$






## Import Curves

Good looking curve data are not very easy to find! I found a few ways to get this data, 

* By doing boundary filtering on black and white (in essense any image with large area of same color could be binarize into black and white) images. `bwboundaries` in `matlab`. 
* Or, we can fetch from raw data of SVG files and any kinds of vector graphics, which contains the coordinates of points on curve. 
  * This is not simple in `matlab`, requires some extra files, like [`loadsvg`](https://ww2.mathworks.cn/matlabcentral/fileexchange/72225-load-svg-into-your-matlab-code), these codes will load the `SVG` files into an cell array containing the coordinates of the geometric object. Add different objects (like different characters will end up in different cells)
  * Note, to prepare `SVG` you may use `Adobe Illustrator` to draw your figure or write the sentence. And then process it to the matlab format. 





