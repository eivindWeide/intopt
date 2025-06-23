# intopt

Program for solving integer [linear programs](https://en.wikipedia.org/wiki/Linear_programming) using the [simplex algorithm](https://en.wikipedia.org/wiki/Simplex_algorithm) written in C.

***
### Compile and run

     gcc intopt.c -lm
     ./a.out < i3

Gives output:

     result is: 192.000000
     Coefficients are: 19.000000*x0 7.000000*x1 0.000000*x2

***
### Input format

The program reads files with the format

    m        n
    
    c_0       c_1      ...  c_n-1
    
    a_0,0     a_0,1    ...  a_0,n-1
    a_1,0     a_1,1    ...  a_1,n-1
                       ...
    a_m-1,0   a_m-1,1  ...  a_m-1,n-1
    
    b_0       b_1      ...  b_m-1

Where *m* is the number of constraints and *n* is the number of variables. First row is the linear function to maximize. Then there is *m* rows of constraints. Last row is *m* boundaries, one for each constraint.

***
### Example

The file **i3**:

    6 3  
    +9 +3 -7
    +5 -8 +2
    +2 +6 +3
    -7 -7 +7
    -1 +0 +0  
    +0 -1 +0
    +0 +0 -1
    +50 +80 -80 +0 +0 +0

is interpreted like 

> Find integers *x0*, *x1*, *x2* to maximize *z* in the function:
>
>>     z = 9x0 + 3x1 - 7x2
>
> While
> 
>>     + 5x0 - 8x1 + 2x2 <= 50
>>     + 2x0 + 6x1 + 3x2 <= 80
>>     - 7x0 - 7x1 + 7x2 <= -80
>>     -x0 <= 0
>>     -x1 <= 0
>>     -x2 <= 0
>
> holds true.
