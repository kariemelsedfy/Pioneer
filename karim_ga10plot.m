%Pioneer
n_0 = 1;
n_1 = 1.46;
n_2 = 2.44;
N = 5;
d_1 = 77.05479452054794;
d_2 = 46.10655737704918;
d_3 = 77.05479452054794;
d_4 = 46.10655737704918;
d_5 = 77.05479452054794;
d_6 = 46.10655737704918;
d_7 = 77.05479452054794;
d_8 = 46.10655737704918;
d_9 = 77.05479452054794;
d_10 = 46.10655737704918;

c = 1;
T_gT = [];
T_gd = [];
    
for lambdas = 370:700 
    B_01 = 1/((2*n_1)) * [(n_1 + n_0) (n_1 - n_0);(n_1 - n_0) (n_1 + n_0)];
    phi_1 = n_1 * ((2 * pi) / lambdas) * d_1;
    P_1 = [exp(-1i*phi_1) 0;0 exp(1i*phi_1)];
    B_12 = 1/((2*n_2)) * [(n_2 + n_1) (n_2 - n_1);(n_2 - n_1) (n_2 + n_1)];
    phi_2 = n_2 * ((2 * pi) / lambdas) * d_2;
    P_2 = [exp(-1i*phi_2) 0;0 exp(1i*phi_2)];
    B_23 = 1/((2*n_1)) * [(n_1 + n_2) (n_1 - n_2);(n_1 - n_2) (n_1 + n_2)];
    phi_3 = n_1 * ((2 * pi) / lambdas) * d_3;
    P_3 = [exp(-1i*phi_3) 0;0 exp(1i*phi_3)];
    B_34 = 1/((2*n_2)) * [(n_2 + n_1) (n_2 - n_1);(n_2 - n_1) (n_2 + n_1)];
    phi_4 = n_2 * ((2 * pi) / lambdas) * d_4;
    P_4 = [exp(-1i*phi_4) 0;0 exp(1i*phi_4)];
    B_45 = 1/((2*n_1)) * [(n_1 + n_2) (n_1 - n_2);(n_1 - n_2) (n_1 + n_2)];
    phi_5 = n_1 * ((2 * pi) / lambdas) * d_5;
    P_5 = [exp(-1i*phi_5) 0;0 exp(1i*phi_5)];
    B_56 = 1/((2*n_2)) * [(n_2 + n_1) (n_2 - n_1);(n_2 - n_1) (n_2 + n_1)];
    phi_6 = n_2 * ((2 * pi) / lambdas) * d_6;
    P_6 = [exp(-1i*phi_6) 0;0 exp(1i*phi_6)];
    B_67 = 1/((2*n_1)) * [(n_1 + n_2) (n_1 - n_2);(n_1 - n_2) (n_1 + n_2)];
    phi_7 = n_1 * ((2 * pi) / lambdas) * d_7;
    P_7 = [exp(-1i*phi_7) 0;0 exp(1i*phi_7)];
    B_78 = 1/((2*n_2)) * [(n_2 + n_1) (n_2 - n_1);(n_2 - n_1) (n_2 + n_1)];
    phi_8 = n_2 * ((2 * pi) / lambdas) * d_8;
    P_8 = [exp(-1i*phi_8) 0;0 exp(1i*phi_8)];
    B_89 = 1/((2*n_1)) * [(n_1 + n_2) (n_1 - n_2);(n_1 - n_2) (n_1 + n_2)];
    phi_9 = n_1 * ((2 * pi) / lambdas) * d_9;
    P_9 = [exp(-1i*phi_9) 0;0 exp(1i*phi_9)];
    B_910 = 1/((2*n_2)) * [(n_2 + n_1) (n_2 - n_1);(n_2 - n_1) (n_2 + n_1)];
    B_inverse = inv(B_89);
    phi_10 = n_2 * ((2 * pi) / lambdas) * d_10;
    P_10 = [exp(-1i*phi_10) 0;0 exp(1i*phi_10)];
    B_100 = 1/((2*n_0)) * [(n_0 + n_2) (n_0 - n_2);(n_0 - n_2) (n_0 + n_2)];
    M =   B_100 *  P_10 *  B_910 * P_9 * B_89 * P_8 * B_78 * P_7 * B_67 * P_6 * B_56 * P_5 * B_45 * P_4 * B_34 * P_3 * B_23 * P_2 * B_12 * P_1 * B_01;
    t = 1/(M(2,2));
    T = abs(t);
    plot(lambdas, T);
    hold on
    T1(c)= T;
    c = c+1;
    %B21_matrix = [AD_2 BC_2 ; BC_2 AD_2];
    %band gap is from 400 to 500
    if (400<lambdas && lambdas<500)
        T_gT = [T_gT;T];
    else
        T_gd = [T_gd;T];
    end
end

T_g = mean(T_gT);
T_d = mean(T_gd);

T_x = T_g/T_d;
T_k = (d_1 + d_2)*N;
hold on
plot(370:700, T1, 'b');
xlabel('Wavelength (nm)')
ylabel('Transmission')