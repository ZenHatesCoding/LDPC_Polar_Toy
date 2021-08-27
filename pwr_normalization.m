function x_out = pwr_normalization(x_in)
   x_out = x_in/sqrt(sum(norm(x_in)^2)/length(x_in));
end