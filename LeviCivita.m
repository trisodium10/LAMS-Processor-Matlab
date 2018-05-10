function epsilon = LeviCivita(i,j,k)
epsilon = -mod((i-j)^2,3)*mod((i-k)^2,3)*mod((j-k)^2,3)*((j-mod(i,3)-1/2)^2-5/4);