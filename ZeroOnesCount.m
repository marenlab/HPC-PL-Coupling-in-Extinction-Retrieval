function [start, len, k1] = ZeroOnesCount(v)
%
% Determine the position and length of each 1-string in v,
% where v is a vector of 1s and 0s
% Derek O'Connor 21 Sep 2011
%
% Based loosely on the WordCount procedure in Kernighan & Plauger [1981], 
% Software Tools in Pascal, ©Bell Labs, Addison-Wesley, page 17.

 n = length(v);
 start = zeros(1,n);            % where each 1-string starts
 len = zeros(1,n);              % length of each 1-string
 k1= 0;                         % count of 1-strings
 inOnes = false;

 for k = 1:n
     if v(k) == 0               % not in 1-string
         inOnes = false;
     elseif ~inOnes             % start of new 1-string
         inOnes = true;
         k1 = k1+1;
         start(k1) = k;
         len(k1) = 1;
     else                       % still in 1-string
         len(k1) = len(k1)+1;
     end

 start(start == 0) = [];
 len(len == 0)     = [];

 end
 
%-------------  ZeroOnesCount ----------------------------
% v = [0 0 0 0 0 0 1 1 1 0 0 0 0 0 1 1 0 0 0 0 1 1 1 1 1 0 0];
% [start,len,k1] = ZeroOnesCount(v);disp([start(1:k1); len(1:k1)]);
%    7    15    21
%    3     2     5