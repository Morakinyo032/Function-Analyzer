# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 19:21:31 2025

@author: HP
"""

import sympy as sp
import numpy as np
from sympy.parsing.sympy_parser import parse_expr

sp.init_printing(use_latex=True)

def is_nested(func: str)-> int:
    '''
    Counts the number of nested absolute value
    symbol that is contained in the function string.

    Parameters
    ----------
    func : str
        String representation of the function.

    Returns
    -------
    int
        Number of nested absolute value symbols.

    '''
    func = func.replace(' ', '')
    stack = []
    length = len(func)
    for i in range(length):
        if func[i] == '|':
            if not stack:
                stack.append(1)
            elif func[i - 1] == '(':
                stack.append(1)
            elif func[i + 1] and func[i + 1]== ')':
                stack.append(-1)
            elif stack[-1] == 1:
                stack.append(-1)
            elif stack[-1] == -1:
                stack.append(1)
    length = len(stack)
    sum = 0
    for i in range(1, length):
        if stack[i] == 1 and stack[i - 1] == 1:
            sum += 1
    if sum == 0:
        return 0
    return sum + 1
        
def parse_abs(func:str)->str:
    '''
    Replaces '|' in the function string to 'abs()'

    Parameters
    ----------
    func : str
        Function String.

    Returns
    -------
    str
        The modified function string.

    '''
    count = is_nested(func)
    if count == 0:
        count = 1
    if '|' in func:
        func = func.replace('|', 'abs(', count=count)
    if '|' in func:
        func = func.replace('|', ')', count=count)
    #print(func)
    return func

class FunctionAnalyser:

    def __init__(self, func_str: str, variable:str, at: float):
        '''
        Initializes an instance of the FunctionAnalyser class

        Parameters
        ----------
        func_str : str
            Function string.
        variable : str
            The variable in the function.
        at : float
            The value of the variable at which the function is be analysed.

        Returns
        -------
        None.

        '''
        self.func_str = func_str
        self.at = at
        self.variable = sp.symbols(variable)
        self.sym_func = self.parse_func_str()
        self.Analysis = f"""
            Let f = {self.func_str}\n
            Funtion f {'is' if self.is_defined() == True else 'is not'}  defined at {self.variable} = {self.at}\n
            The value of f at {self.variable} = {self.at} is {self.solve_at()}\n
            The left hand side limit of f at {self.variable} = {self.at} is {self.left_limit()}\n
            The right hand side limit of f at {self.variable} = {self.at} is {self.right_limit()}\n
            Funtion f {'is' if self.is_continuous() == True else 'is not'}  continuous at {self.variable} = {self.at}\n
            The left hand side derivative of f at {self.variable} = {self.at} is {self.left_deriv()}\n
            The right hand side derivative of f at {self.variable} = {self.at} is {self.right_deriv()}\n
            Funtion f {'is' if self.is_differentiable() == True else 'is not'}  differentiable at {self.variable} = {self.at}\n
            """
        
    def parse_func_str(self, func=None)-> sp.EX:
        '''
        Parses the function string to sympy experession.

        Parameters
        ----------
        func : str, optional
            The default is None.

        Returns
        -------
        sp.Ex
            Sympy representation of the function string.

        '''
        try:       
            if func == None:
                func = self.func_str
            while '|' in func:
                func = parse_abs(func)
            return sp.sympify(func)
        
        except Exception:
            print('Cannot parse the function string.\
                  Ensure all nested "|" are properly enclosed in brackects')
            
            
    def is_defined(self, func=None, value=None)-> tuple:
        '''
        Check if a function is defined at the specified point.

        Parameters
        ----------
        func : sp.EX, optional
            The default is None.
        value : float, optional
           The default is None.

        Returns
        -------
        tuple
            Tuple containing true or false value and the value of the function at that point.

        '''
        if func == None:
            func = self.sym_func
        if value == None:
            value = self.at
        result = func.subs(self.variable, value).evalf()
        #print(result)
        if result.is_finite:
            return True
        return False
    
    def solve_at(self, func=None, variable=None, at=None)-> float:
        '''
        Solves the given function at a specified value

        Parameters
        ----------
        func : sp.Ex , optional
            DESCRIPTION. The default is None.
        variable : str, optional
            DESCRIPTION. The default is None.
        at : float, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        float
            The value of the function at the given point.

        '''
        if variable == None:
            variable = self.variable
        if func == None:
            func = self.sym_func
        if at == None:
            at = self.at
        value = func.subs(variable, at).evalf()
        #print(value)
        if value.is_finite:
            return value
        return np.inf
        
    def left_limit(self, func=None, variable=None, at=None):
        '''
        Claculate the left-hand side limit of the given function.

        Parameters
        ----------
        func : sp.Ex, optional
            DESCRIPTION. The default is None.
        variable : str, optional
            DESCRIPTION. The default is None.
        at : float, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        float
            The limit of the function from the left.

        '''
        if variable == None:
            variable = self.variable
        if func == None:
            func = self.sym_func
        if at == None:
            at = self.at
        llimit = sp.limit(func, variable, at, '-').evalf()
        if llimit.is_finite:
            return llimit
        return np.inf
        
    def right_limit(self, func=None, variable=None, at=None):
        '''
        Claculate the right-hand side limit of the given function.

        Parameters
        ----------
        func : sp.Ex, optional
            DESCRIPTION. The default is None.
        variable : str, optional
            DESCRIPTION. The default is None.
        at : float, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        float
            The limit of the function from the left.

        '''
        if variable == None:
            variable = self.variable
        if func == None:
            func = self.sym_func
        if at == None:
            at = self.at
        rlimit = sp.limit(func, variable, at, '+').evalf()
        if rlimit.is_finite:
            return rlimit
        return np.inf
    
    def f_prime_at(self, func=None, variable=None, at=None):
        '''
        Calculate the value of the derivate of the given function at the given point.

        Parameters
        ----------
        func : sp.Ex, optional
            DESCRIPTION. The default is None.
        variable : str, optional
            DESCRIPTION. The default is None.
        at : float, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        float
            The derivative of the function st the given point.

        '''
        if variable == None:
            variable = self.variable
        if func == None:
            func = self.sym_func
        if at == None:
            at = self.at
        derivative = sp.diff(func, variable)
        return derivative

    def is_continuous(self, func=None, variable=None, at=None):
        '''
        Checks if a function is continuous at a given point.

        Parameters
        ----------
        func : sp.Ex, optional
            DESCRIPTION. The default is None.
        variable : str, optional
            DESCRIPTION. The default is None.
        at : float, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        bool
            True if the function is continuous, else, false.

        '''
        if variable == None:
            variable = self.variable
        if func == None:
            func = self.sym_func
        if at == None:
            at = self.at
        define = self.is_defined()
        if define and (self.left_limit(func, variable, at) == self.right_limit(func, variable, at) == self.solve_at(func, variable, at)):
            return True
        return False
    
    def left_deriv(self, func=None, variable=None, at=None):
        '''
        Calculate the left-hand side derivative of the given function at a given point.

        Parameters
        ----------
        func : sp.Ex, optional
            DESCRIPTION. The default is None.
        variable : str, optional
            DESCRIPTION. The default is None.
        at : float, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        float
            The derivative of the function from the left.

        '''
        if variable == None:
            variable = self.variable
        if func == None:
            func = self.sym_func
        if at == None:
            at = self.at
        func_at_point = func.subs(variable, at)
        func = sp.sympify((func - func_at_point)/(variable - at))
        lderiv = sp.limit(func, variable, at, '-').evalf()
        if lderiv.is_finite:
            return lderiv
        return np.inf
    
    def right_deriv(self, func=None, variable=None, at=None):
        '''
        Calculate the right-hand side derivative of the given function at a given point.

        Parameters
        ----------
        func : sp.Ex, optional
            DESCRIPTION. The default is None.
        variable : str, optional
            DESCRIPTION. The default is None.
        at : float, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        float
            The derivative of the function from the right.

        '''
        if variable == None:
            variable = self.variable
        if func == None:
            func = self.sym_func
        if at == None:
            at = self.at
        func_at_point = func.subs(variable, at)
        func = sp.sympify((func - func_at_point)/(variable - at))
        rderiv = sp.limit(func, variable, at, '+').evalf()
        if rderiv.is_finite:
            return rderiv
        return np.inf

    def is_differentiable(self, func=None, variable=None, at=None)-> bool:
        '''
        Checks if a function is differentiable at a given point.

        Parameters
        ----------
        func : sp.Ex, optional
            DESCRIPTION. The default is None.
        variable : str, optional
            DESCRIPTION. The default is None.
        at : float, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        bool
            True if the function is differentiable, else, false.

        '''
        
        if variable == None:
            variable = self.variable
        if func == None:
            func = self.sym_func
        if at == None:
            at = self.at        
        if self.is_continuous(func, variable, at) and (self.left_deriv(func, variable, at) == self.right_deriv(func, variable, at)) and self.left_deriv(func, variable, at) != np.inf:
            return True
        return False

function = 'exp(-x^2)'
analyse = FunctionAnalyser(function,'x', 100)
print(analyse.Analysis)
expr = analyse.f_prime_at()
sp.pprint(expr)
