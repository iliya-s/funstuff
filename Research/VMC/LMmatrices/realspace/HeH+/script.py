import numpy as np

d = 10
#x = np.zeros(d)
x = np.random.rand(d)
print("sin(\sum_i^d x_i)")
lhs = np.sin(np.sum(x))
print(lhs)
print("\t")

print("\sum_j^d \prod_k,k!=j^d sin(x_k + alpha_k - alpha_j) / sin(alpha_k - alpha_j)")
rhs = 0.0
alpha = np.random.rand(d)
for j in range(len(x)):
    val = np.sin(x[j])
    for k in range(len(x)):
        if k != j:
            val *= np.sin(x[k] + alpha[k] - alpha[j]) / np.sin(alpha[k] - alpha[j])
    rhs += val;

print(rhs)
