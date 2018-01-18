# cube-gpu-raytracing
GPU Ray tracing for the world which only exists cube

## Result

### Render

256 Samples : 

![256 samples](out%20(14).bmp)

### Performance

#### CPU Threading

| Samples\Threading | 1 | 2 | 4 | 8 | 16 |
|:--:|-----:|-----:|-----:|-----:|-----:|
| 1  | 1.42 | 0.82 | 0.44 | 0.25 | 0.21 |
| 4  | 5.28 | 2.74 | 2.08 | 1.11 | 0.62 |
| 9  | 11.95| 6.27 | 3.45 | 1.81 | 1.29 |
|16  | 21.00|11.41 | 5.96 | 4.28 | 2.22 |
|25  | 32.46|18.50 | 9.64 | 4.72 | 3.44 |

#### CPU vs GPU

| Samples\Type | CPU(1 Thread) |  CPU(16 Threads) | GPU |
|:---:|-----:|-----:|-----:|
| 1   |1.42  | 0.21 | 0.23 |
| 4   |5.28  | 0.61 | 0.85 |
| 16  |21.00 | 2.26 |3.24 |
|256  |348.10| 23.86|51.21|
