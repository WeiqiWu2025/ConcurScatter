function [H_LS] = funcLSEstimator(x,y)
% 输入
    % x: 传输信号
    % y: 接收信号
% 输出
    % H_LS: LS信道估计
    
    H_LS = inv(x'*x)*x'*y;
end


