Bez(t, p1, p2, a, b) = (1-t)^3 * p1 + 3 * t * (1-t)^2 * a + 3 * t^2 * (1-t) * b + t^3 * p2; 

function get_BezMat(N)
    Bez_mat = zeros(N,N);
    for n in 2:(N-1)
        if(n>1 && n<N)
            Bez_mat[n,n] = 4; Bez_mat[n,n+1] = 1; Bez_mat[n,n-1] = 1;
        end
    end
    Bez_mat[1,1] = 2; Bez_mat[1,2] = 1;
    Bez_mat[N,N-1] = 2; Bez_mat[N,N] = 7;
    return Bez_mat
end;
function get_BezVec(N, x)
    Bez_vec = zeros(N);
    for n in 2:(N-1)
        Bez_vec[n] = 2*(2*x[n] + x[n+1])
    end
    Bez_vec[1] = x[1] + 2*x[2]
    Bez_vec[N] = 8*x[N] + x[N+1]
    return Bez_vec
end;
function get_a2b_Bez(x,a)
    N = length(a)
    b = zeros(N)
    for n in 1:(N-1)
        b[n] = 2*x[n+1] - a[n+1]
    end
    b[N] = 0.5*(a[N] + x[N+1])
    return b
end;

# Should be length(p_set) = N+1
function get_a_b_Bezier(N, p_set)
    Bez_mat = get_BezMat(N)
    Bez_vec = get_BezVec(N, p_set);
    
    #Ax=b
    #x = cholesky(A) \ b # if A is Hermetian
    a_vec = Bez_mat \ Bez_vec;
    b_vec = get_a2b_Bez(N, p_set, a_vec);
    return (a_vec, b_vec)
end; 
function get_integrated_Bezier_2nd()
    a = 1.0/7
    b = 1.0/42
    c = 1.0/105
    d = 1.0/140
    return B2 = Matrix([[a, b, c, d] [b, c, d, c] [c, d, c, b] [d, c, b, a]])
end;


function get_Bezier_vec_Mat_for_all_sites_time(L, x1_traject, x2_traject)
    #--Interpolation points for single site freq. for each time --#
    n_traject = size(x1_traject,1)
    N = n_traject-1
    point_set_x1 = zeros(L, n_traject)
    point_set_x1_a = zeros(L, N)
    point_set_x1_b = zeros(L, N)
    for i in 1:L
        point_set_x1[i,:] = [ x1_traject[n][i] for n in 1:n_traject];
        size(point_set_x1), size(point_set_x1[i,:]), N
        (a_vec_i, b_vec_i) = get_a_b_Bezier(N, point_set_x1[i,:])
        #@show size(a_vec_i), size(b_vec_i)
        point_set_x1_a[i,:] = copy(a_vec_i)
        point_set_x1_b[i,:] = copy(b_vec_i)
    end
    #--------------------------#
    #-- Interpolation points for corr. for each time --#
    point_set_x2 = zeros(L,L, n_traject)
    point_set_x2_a = zeros(L,L, N)
    point_set_x2_b = zeros(L,L, N)

    for i in 1:L
        for j in i:L
            point_set_x2[i,j,:] = [ x2_traject[n][i,j] for n in 1:n_traject];
            (a_vec_i, b_vec_i) = get_a_b_Bezier(N, point_set_x2[i,j,:])
            point_set_x2_a[i,j,:] = copy(a_vec_i)
            point_set_x2_a[j,i,:] = copy(a_vec_i)
            point_set_x2_b[j,i,:] = copy(b_vec_i)
            point_set_x2_b[i,j,:] = copy(b_vec_i)

        end
    end
    return (point_set_x1_a, point_set_x1_b, point_set_x2_a, point_set_x2_b)
end;
