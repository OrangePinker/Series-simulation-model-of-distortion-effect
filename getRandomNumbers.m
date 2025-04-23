function [num1, num2] = getRandomNumbers()
    nums = randperm(512, 2);
    num1 = nums(1);
    num2 = nums(2);
end