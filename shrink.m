%Shrink function as specified in Kilmer (2011)
function y = shrink(x,epsilon)
y = 0;
if x > epsilon
    y = x-epsilon;
elseif x < -epsilon
    y = x+epsilon;
end

end