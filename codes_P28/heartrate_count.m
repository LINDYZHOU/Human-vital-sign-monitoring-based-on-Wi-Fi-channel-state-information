h_data = importdata("ax40.txt");
h_data_count = h_data.data(10:end-5,1);
fprintf('mean: %4.3f, median: %4.3f, mode: %4.3f.\n', mean(h_data_count), ...
    median(h_data_count),mode(h_data_count));
