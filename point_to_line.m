function d = point_to_line(pt, fit)  % y = fit_1 *x + fit_2; pt=  [pt_x, pt_y]
      fit_1 = fit(1);
      fit_2 = fit(2);
      numr = abs(fit_1*pt(1)+fit_2-pt(2));  % numerator
      denr = sqrt(fit_1^2+1);  % denominator;
      d = numr/denr; %distance from point pt to line fit_1 
end