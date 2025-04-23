% 定义列号转Excel列字母的函数
function column_letter = get_column_letter(col_num)
    letters = '';
    while col_num > 0
        remainder = mod(col_num - 1, 26);
        letters = [char(65 + remainder), letters]; % 65对应'A'
        col_num = floor((col_num - 1) / 26);
    end
    column_letter = letters;
end