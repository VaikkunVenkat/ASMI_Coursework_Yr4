function y = sinegen(fsamp,fsig,nsamp)
    tsamp = 1/fsamp;
    t = 0 : tsamp : (nsamp-1)*tsamp;
    y = sin(2*pi*fsig*t);
end

