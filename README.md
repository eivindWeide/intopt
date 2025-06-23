# intopt

[Simplex](https://en.wikipedia.org/wiki/Simplex_algorithm) program for solving integer [linear programs](https://en.wikipedia.org/wiki/Linear_programming)

#### Compile and run

> gcc intopt.c -lm

> ./a.out < i3

Gives output

> result is: 192.000000

> Coefficients are: 19.000000\*x0 7.000000\*x1 0.000000\*x2

#### Input format

The program reads files with the format

> m  n
> c0  c1  ...  xn-1
> a0,0  a0,1  ...  a0,n-1
> a1,0  a1,1  ... a1,n-1
>   ...
> am-1,0  am-1,1  ...  am-1,n-1
> b0  b1  ...  bm-1

Where *m* is the number of constraints and *n* is the number of decision variables. First row is the linear function to maximize. Then there is m rows of constraints. Last row is the bound for each constraint.

The file `i3`

> 6 3
> +9 +3 -7
> +5 -8 +2
> +2 +6 +3
> -7 -7 +7
> -1 +0 +0 
> +0 -1 +0
> +0 +0 -1
> +50 +80 -80 +0 +0 +0

is interprited as 

> Maximize 9x0 + 3x1 - 7x2
> While
>   5x0 - 8x1 + 2x2 <= 50
>   2x0 + 6x1 + 3x2 <= 80
>   - 7x0 - 7x1 + 7x2 <= -80
>   - 1x0 + 0x1 + 0x1 <= 0  (x >= 0)  
>   0x0 - 1x1 + 0x2 <= 0
>   0x0 + 0x1 - 1x2 <= 0
> holds true
