function coord = coordi_CT(x_in, z_in, img, img_start)
    coord(1) = (x_in + img_start(1))*img.spacing(1)+img.origin(1);
    coord(2) = img.origin(3);
    coord(3) = -(z_in+img_start(2))*img.spacing(2)+img.origin(2);
end