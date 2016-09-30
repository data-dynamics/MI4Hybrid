function output = NUM2STR(i)

if i >=0
    output = num2str(i);
else
    output = ['(- ' num2str(abs(i)) ')'];
end