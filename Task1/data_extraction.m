function [speed_x, speed_y, pressure,position_x,position_y] = data_extraction(A1)
    
speed_x = A1(1,:);
    speed_y = A1(2,:);
    pressure = A1(3,:);
    position_x = cumsum(speed_x);
    position_y = cumsum(speed_y);
    
end