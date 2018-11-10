import tensorflow as tf
import numpy as np
import conjgrad as c

#Analytical result
print("Analytic result")
xa = 1.
ya = 2.
fa = xa * np.sin(ya) + (xa / ya)
dfa_dxa = np.sin(ya) + 1 / ya
dfa_dya = xa * np.cos(ya) - xa / (ya*ya)
print(fa)
print(dfa_dxa, dfa_dya)


#Tensorflow example#
print("Tensorflow example")
#define tensorflow variables
x = tf.placeholder(tf.float64)
y = tf.placeholder(tf.float64)

#define function
f = x * tf.sin(y) + (x / y)

#create gradient object
grad = tf.gradients(f, [x, y])

#create session
s = tf.Session()
val = s.run([f, grad], {x: [1], y: [2]} )
s.close

print(np.array(val[0]))
print(np.array(val[1]))

#Tensorflow matrix example#
print("Tensorflow matrix example")
#A*x = b, calculating derivative dx_dA(i,j) = sum_l dx(l)_dA(i,j)..... I think
    #returns sum(dy/dx) for each independent variable, x
#Build tensorflow objects
A = tf.placeholder(dtype = tf.float32, shape=(3,3)) #A matrix of placeholder vals
b = tf.constant([[1.0],[2.0],[3.0]])    #b vector
#x = A^-1 * b
x = tf.matmul(tf.matrix_inverse(A),b)

grad = tf.gradients(x, A)

s = tf.Session()
val = s.run([x, grad], {A : [[5.0,1.0,1.0],[1.0,10.0,1.0],[1.0,1.0,15.0]]})
s.close

print(np.array(val[0]))
print(np.array(val[1]))

#Try to do matrix derivative with iterative function
print("My Tensorflow matrix example")
A = tf.placeholder(dtype=tf.float64, shape = (3,3))
b = tf.constant([[1.0],[2.0],[3.0]],dtype=tf.float64)

x0 = tf.random_uniform((3,1),dtype=tf.float64)
r0 = b - tf.matmul(A,x0)
p0 = r0
dotr0 = tf.matmul(r0, r0, transpose_a = True)
for i in range(30):
    Ap = tf.matmul(A,p0)
    pAp = tf.matmul(p0, Ap, transpose_a = True)
    a = dotr0/pAp
    x1 = x0 + tf.scalar_mul(tf.reshape(a,[]),p0)
    r1 = r0 - tf.scalar_mul(tf.reshape(a,[]),Ap)
    dotr1 = tf.matmul(r1, r1, transpose_a = True)
    b = dotr1/dotr0
    p1 = r1 + tf.scalar_mul(tf.reshape(b,[]),p0)
    x0 = x1
    r0 = r1
    p0 = p1
    dotr0 = dotr1

grad = tf.gradients(x1, A)

s = tf.Session()
val = s.run([x1, grad], {A : [[5.0,1.0,1.0],[1.0,10.0,1.0],[1.0,1.0,15.0]]})
s.close

print(np.array(val[0]))
print(np.array(val[1]))

