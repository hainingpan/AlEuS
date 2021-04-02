% f()
% 
% function f()
% disp('2')
% 
% 
% end
% 
% function h()
% disp('3')
% end

clear all
% for i = 1:50
    x = rand(2000,2000);
    tic;

    [V,D]=eig(x);
    toc;

% end
