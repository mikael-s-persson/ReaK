function count_out = add_predstats(fName,outName)

cur_s = load(fName);
[cur_n, cur_m] = size(cur_s);

prev_s = load(outName);
[prev_n, prev_m] = size(prev_s);

if(prev_m != cur_m+1)
    return
end

if(prev_n < cur_n)
    prev_s = [prev_s; zeros(cur_n - prev_n, prev_m)];
end

for I=1:cur_n
    prev_s(I,1) = cur_s(I,1);
    cur_s_fin = isfinite(cur_s(I,3:cur_m));
    if(sum(cur_s_fin) == cur_m - 2) % if they are all finite numbers
        prev_s(I,3:prev_m) = (prev_s(I,2) * prev_s(I,3:prev_m) + cur_s(I,2:cur_m)) / (prev_s(I,2) + 1);
        prev_s(I,2) = prev_s(I,2) + 1;
    end
end

save(outName, 'prev_s', '-ascii');

end