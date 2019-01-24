︠91cf517a-d6e4-42e4-8308-0faa9a8d3f76s︠
sage_server.MAX_STDOUT_SIZE=sage_server.MAX_STDOUT_SIZE*10
import itertools
q=128
o=2
v=2*o
n=v+o
J=AffineGroup(n, GF(q))
K.<a>=GF(q)
K2.<a>=GF(2)
Pp=PolynomialRing(K,'x',n)
P2=PolynomialRing(K2,'x',n)
x=Pp.gens()
x_vec=vector(x)
tp1 = matrix(x)
m=13
#----------------------------------------------------Public key generation
A00=[0 for i in range(o)]
A01=[0 for i in range(o)]
A02=[0 for i in range(o)]
A03=[0 for i in range(o)]

B0=[0 for i in range(o)]
B=[[] for i in range(o)]
A=[0 for i in range(o)]
N=[0 for i in range(o)]

D=[0 for i in range(o)]
D[0]=random_matrix(K2,v,v) #--- vvpart

L=[0 for i in range(o)]
L[0]=random_matrix(K2,v,o) #---vo part

I=[0 for i in range(o)]
I[0]=random_matrix(K2,o,v) #---ov part

M=[0 for i in range(o)]
M[0]=matrix(K2,o,o)


for i in range(o):  #------D[i] is the vv part
    D[i] = random_matrix(K2,v,v)
    #parent(D[i])

for i in range(o):  #-------L[i] is the vo part
    L[i] = random_matrix(K2,v,o)

for i in range(o): #---------I[i] is the ov part
    I[i] = random_matrix(K2,o,v)

for i in range(o): #--------M[i] is the oo part
    M[i] = matrix(K2,o,o)

for i in range(o):
    A00[i]=D[i]
    A01[i]=L[i]
    A02[i]=I[i]
    A03[i]=M[i]

for i in range(o):
    N[i] = matrix.block([[A00[i],A01[i]],[A02[i],A03[i]]])

F=[0 for i in range(o)]
print "Generating Oil-Vinegar matrices..."
for i in range(o):
    F[i]=N[i]#-------------------central quadratic matrix
    print F[i]
    print""
    tp1*F[i]*tp1.transpose()
    print ""

while true:
        T=random_matrix(K2,n)
        if T.is_invertible():
            break
Tt=T.transpose()

P=[0 for i in range(o)]
print "Generating public key polynomials... "
for i in range(o):
    P[i]=Tt*(F[i].change_ring(K))*T
    #parent(P[i])
    print tp1*P[i]*tp1.transpose()
print ""
#-------------------------------------------------------------------------------------Signature generation
coremap = [tp1*F[i]*tp1.transpose() for i in range(o)]
print "Core map F..."
for i in range(o):
    F[i]=N[i]
    print tp1*F[i]*tp1.transpose()
Ry=PolynomialRing(K,n-v,['x%s'%p for p in[v..n-1]])
images = [K.random_element() for i in range(v)]+list(Ry.gens())
phi=Pp.hom(images,Ry)
List=[coremap[j] for j in  range(o)]
coremap_subs = [phi(f[0]) for f in List]
print ""
def hash_function(m):
    m_minus_one = m-1
    list = []
    for i in range(o):
        list.append(ZZ.random_element(m))
    return list

hashresult=hash_function(m)
print("Hash result:")
print(hashresult)
print ""

s = [coremap_subs[i]-hashresult[i] for i in range(o)]
h = Ry.ideal(s)

def give_result(ring,coremap):
    h = ring.ideal(coremap)
    print "Solving for oil variables..."
    result = h.variety()[0]
    print(result)
    print("")
    for i in [v..n-1]:
        images[i] = result['x'+str(i)]
    print "Pre-image of the hash value under the core map F..."
    print images
    return images

pre_image = give_result(Ry,s)
preimage = vector(pre_image)

signature = T.inverse()*preimage
print ""
print "Computing signature..."
print signature
print ""
#----------------------------------------------------------------------Signature verification
publickey = [tp1*P[i]*tp1.transpose() for i in range(o)]
#Rz=PolynomialRing(K,n-v,['x%s'%p for p in [v..n-1]])
#Rz=PolynomialRing(K,n,['x%s'%p for p in[0..n-1]])
images = [signature[i] for i in range(n)]
phi=Pp.hom(images,K)
myList=[publickey[j] for j in  range(o)] 
publickey_subs = [phi(f[0]) for f in myList]
print "Signature validation..."
publickey_subs == hashresult
all(itertools.imap(lambda x: x in publickey_subs, hashresult))





