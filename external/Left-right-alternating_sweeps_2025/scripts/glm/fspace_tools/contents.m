% This folder contains a set of tools for working with function spaces,
% particularly for using function spaces as a means of decomposing
% variables (functional decomposition).
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESSENTIAL COMPONENTS OF THIS TOOLKIT:
% 
% BASISFUNCTION - class respresenting basis functions, the basis vectors of
% a function space.
%
% FUNCTIONSPACE - a class representing a collection of BasisFunction
% objects. This class provides useful methods for performing and
% manipulating functional decompositions of variables
%
% DRFUNCTIONSPACE - an extension of FuntionSpace which dimension-reduces
% the function output
%
% FDECOMP - a class representing a functional decomposition of a variable
%
% CREATEFCNSPACE - a convenience function for creating a FunctionSpace
% using preset options