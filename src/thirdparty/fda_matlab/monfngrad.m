function [f, grad] = monfngrad(x, Wfd)
%MONFNGRAD evaluates a monotone function of the form
%  h(x) = (D^{-1} exp Wfd)(x)
%  where  D^{-1} means taking the indefinite integral,
%  and its gradient

%  Last modified 10 December 2000

f    =   monfn(x, Wfd);
grad = mongrad(x, Wfd);
