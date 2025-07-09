%% colours and line formatting
colours = [ 0.0, 0.7, 0.8;
           0.12, 0.34, 0.57;
           0.91, 0.15, 0.76;
           0.31, 0.12, 0.77;
           0.93, 0.13, 0.65;
           0.55, 0.51, 0.87;
           0.61, 0.78, 0.79;
           0.01, 0.31, 0.39;
           0.71, 0.25, 0.81;
           0.83, 0.69, 0.44;
           0.06, 0.40, 0.74;
           0.18, 0.18, 0.53;
           0.34, 0.72, 0.53;
           0.94, 0.38, 0.64;
           0.70, 0.15, 0.88;
           0.60, 0.67, 0.09;
           0.91, 0.29, 0.31;
           0.80, 0.86, 0.31;
           0.19, 0.93, 0.42;
           0.95, 0.79, 0.21;
           0.14, 0.41, 0.05
         ];

% Define line styles and markers for plotting
linestyles = {'-', '--', ':', '-.'};
markers = {'o', 's', 'd', '^'};

%% Starting

baseGraph5GNR = 'NR_1_5_352'; 
coderate = [1/3 1/2 3/5 4/5];  
eb_no_dbvec = 0:0.5:10;
[B,Hfull,z] = nrldpc_Hmatrix(baseGraph5GNR,352); % Convert the base H matrix to binary H matrix with z=352

% Simulation parameters
nsim = 1;  
max_it = 50; 
iterations = 1:1:max_it;

% Arrays to store simulation results for all code rates
success_prob_vs_ebno = zeros(length(coderate), length(eb_no_dbvec));
ber_vs_ebno = zeros(length(coderate), length(eb_no_dbvec));

figure(100); % Figure for Success Probability vs Eb/No plots
hold on;

%% For different coderates and Monte Carlo

cr_idx = 1;
for cr = coderate

    % Create a new figure for success probability vs iteration for this code rate
    figure(cr_idx);
    hold on;

    % Performing rate matching
    [mb,nb] = size(B); % mb=46, nb=68 for NR_1_5_352
    kb = nb - mb; % Number of information columns: 68-46=22
    kNumInfoBits = kb * z; % Number of information bits: 22*352=7744

    k_pc = kb-2; % Systematic bits excluding 2 punctured bits: 22-2=20
    nbRM = ceil(k_pc/cr); % Number of coded bits per z, adjusted for code rate
    nBlockLength = nbRM * z; % Number of encoded bits

    H = Hfull(:,1:nBlockLength);
    nChecksNotPunctured = mb*z - nb*z + nBlockLength;
    if nChecksNotPunctured > size(Hfull,1)
        nChecksNotPunctured = size(Hfull,1); % Ensure no row index violation
    end
    H = H(1:nChecksNotPunctured,:); % Binary H matrix after puncturing

    % Rate matching done
    [row,col] = size(H);
    L = zeros(size(H));
    k = col - row;

    % Performing soft decoding
    cn_to_vn_map = cn_vn(H); % Shows ith cn connected to which all vns
    vn_to_cn_map = vn_cn(H); % Shows ith vn connected to which all cns

    d_iter = 1;
    decoding_error = zeros(1,length(eb_no_dbvec));
    bit_error = zeros(1,length(eb_no_dbvec));

    eb_no_idx = 1;
    for eb_no_db = eb_no_dbvec
        eb_no = 10^(eb_no_db/10);
        sigma = sqrt(1/(2*cr*eb_no));
        success = 0;
        error1 = 0;  % Total bit errors
        total_bits = 0;  % Total bits transmitted
        itr_success = nsim.*ones(1,max_it);
        vn_sum_vec = zeros(1,col);

        for sim=1:nsim

            org_msg = randi([0 1],[k 1]); % Generate information bit vector
            encoded_msg = nrldpc_encode(B,z,org_msg'); % Encode using 5G NR LDPC base matrix
            encoded_msg = encoded_msg(1:nBlockLength);

            n = length(encoded_msg);
            % Performing BPSK modulation
            bpsk_msg = 1 - 2.*encoded_msg;
            % Generating noise
            noise = sigma * randn(1,n);

            received_bpsk = bpsk_msg + noise;
            % Compute channel LLRs for sum-product decoding
            llr = (2 / sigma^2) * received_bpsk;
            % Changing message back to bits for received msg only
            received_bits = (received_bpsk<0);
            prev_msg = received_bits;
            c_hat = zeros(1,col);

            for it = 1:max_it

                % Message from VN to CN
                if(it==1)
                    for i=1:col
                        for j=vn_to_cn_map{i,1}
                            L(j,i) = llr(1,i);
                        end
                    end
                else
                    for i = 1:col
                        for j=vn_to_cn_map{i,1}
                            L(j,i) = vn_sum_vec(1,i) - L(j,i);
                        end
                    end
                end

                % Message from CN to VN using sum-product algorithm
                for i=1:row
                    for j=cn_to_vn_map{i,1}
                        prod_term = 1;
                        for jj=cn_to_vn_map{i,1}
                            if jj~=j
                                prod_term = prod_term * tanh(L(i,jj)/2);
                            end
                        end
                        L(i,j) = 2 * atanh(prod_term);
                    end
                end

                % Finding sum of values received by each VN
                for i = 1:col
                    sum1 = llr(1,i);
                    temp = L(:,i);
                    sum1 = sum1 + sum(temp);
                    vn_sum_vec(1,i)=sum1;
                end

                c_hat = (vn_sum_vec<0);

                if(sum(xor(c_hat(1:k),org_msg'))==0)
                    success = success+1;
                    break;
                else
                    itr_success(1,it)=itr_success(1,it)-1;
                end

                % Calculate bit errors for BER calculation
                bit_errors_this_frame = sum(c_hat ~= encoded_msg);
                error1 = error1 + bit_errors_this_frame;
                total_bits = total_bits + col;

                if(sum(xor(prev_msg,c_hat))==0)
                    for tmp_itr=it+1:max_it
                        itr_success(1,tmp_itr)=itr_success(1,tmp_itr)-1;
                    end
                    break;
                end
                prev_msg = c_hat;

            end

        end

        figure(cr_idx);
        bx = gca;
        bx.YScale = 'log';
        bx.YLim   = [1e-3 1];
        bx.YTick  = 10.^(-5:0);
        plot(iterations, itr_success./nsim, 'Color', colours(eb_no_idx,:), 'DisplayName', [num2str(eb_no_db) ' dB']);

        decoding_error(1, eb_no_idx) = (nsim-success)/nsim;
        bit_error(1, eb_no_idx) = error1/total_bits;

        % Store results in our arrays for later plotting
        success_prob_vs_ebno(cr_idx, eb_no_idx) = success/nsim;  % Success probability
        ber_vs_ebno(cr_idx, eb_no_idx) = error1/total_bits;      % BER

        eb_no_idx = eb_no_idx + 1;
    end

    % Finalize the Success Probability vs Iteration plot for this code rate
    figure(cr_idx);
    xlabel("Iteration number");
    ylabel('Success Probability at each iteration');
    title(['Success Probability vs. Iteration for Soft Decoding, Code Rate = ', num2str(cr)]);
    grid on;
    legend('show');
    hold off;

    % Plot Success Probability vs Eb/No for this code rate in the main figure
    figure(100);
    plot(eb_no_dbvec, success_prob_vs_ebno(cr_idx, :), 'LineWidth', 2, 'DisplayName', ['Coderate = ', num2str(cr)]);
    
    cr_idx = cr_idx + 1;
end

% Finalize Success Probability vs Eb/No plot
figure(100);
xlabel("Eb/No (dB)");
ylabel("Success Probability");
title("Success Probability vs Eb/No for Different Code Rates");
legend('show');
grid on;
hold off;

%% Normal approximation function
function [P_err] = normal_approximation(n, k, eb_no_db, rate)
   
    % Convert Eb/N0 from dB to linear
    eb_no = 10.^(eb_no_db/10);
 
    snr = eb_no * rate;
    
    % Channel Capacity for BPSK (binary-input AWGN)
    C = log2(1 + snr);
    
    % Channel Dispersion for BPSK (binary-input AWGN)
    V = (1 - 1 ./ (1 + snr).^2) * (log2(exp(1)))^2;
    
    % Normal approximation bound
    epsilon = qfunc( sqrt(n ./ V) .* (C - rate));
    
    % Return error probability
    P_err = epsilon;
end

% Create a dedicated figure for normal approximation comparison
figure(200);
hold on;

for i = 1:length(coderate)
    % Plot actual BER
    semilogy(eb_no_dbvec, ber_vs_ebno(i,:), '-', 'LineWidth', 2, 'Marker', markers{mod(i-1,4)+1}, 'MarkerSize', 6, ...
        'DisplayName', ['Actual BER, R = ', num2str(coderate(i))]);
    
    % Calculate normal approximation for each Eb/N0 point
    n = nBlockLength_vec(i);  % Use the block length from the simulation
    k = round(coderate(i) * n);  % match actual rate!

    normal_approx_ber = zeros(1, length(eb_no_dbvec));
    
    for j = 1:length(eb_no_dbvec)
        normal_approx_ber(j) = normal_approximation(n, k, eb_no_dbvec(j), coderate(i));
    end
    
    % Plot normal approximation
    semilogy(eb_no_dbvec, normal_approx_ber, '--', 'Color', colours(i,:), 'LineWidth', 1.5, ...
        'DisplayName', ['Normal Approx, R = ', num2str(coderate(i))]);
    
    % Calculate and plot Shannon bound for this code rate
    shannon_bound_ebno_linear = (2^coderate(i) - 1)/(coderate(i));
    shannon_bound_ebno_db = 10*log10(shannon_bound_ebno_linear);
    
    % Plot vertical line at Shannon bound
    xline(shannon_bound_ebno_db, ':', ['Shannon limit, R = ', num2str(coderate(i))], 'LineWidth', 1.5);
end

ax = gca;
ax.YScale = 'log';                   % force log‑scale y
ax.YLim   = [1e-5 1];                
ax.YTick  = 10.^(-5:0);            
xlabel("Eb/No (dB)");
ylabel("Bit Error Rate (BER)");
title("BER Performance with Normal Approximation and Shannon Limits of Sum-Product Algorithm");
legend('show');
ylim([1e-5 1]);
grid on;
hold off;

%% Helper functions

function [B,H,z] = nrldpc_Hmatrix(BG,z)
    % Load the base matrix NR_1_5_352
    % Note: The actual loading of the matrix file depends on your file structure
    % Here, we assume NR_1_5_352 is available as a variable or file
    % Replace the following line with actual loading mechanism if needed
    load(sprintf('%s.txt',BG),BG); % This line assumes a text file exists
    B = eval(BG); % Use the loaded base matrix
    [mb,nb] = size(B);
    H = zeros(mb*z,nb*z); 
    Iz = eye(z); I0 = zeros(z);
    for kk = 1:mb
        tmpvecR = (kk-1)*z+(1:z);
        for kk1 = 1:nb
            tmpvecC = (kk1-1)*z+(1:z);
            if B(kk,kk1) == -1
                H(tmpvecR,tmpvecC) = I0;
            else
                H(tmpvecR,tmpvecC) = circshift(Iz,-B(kk,kk1));
            end
        end
    end
end

function out=cn_vn(H)
    [row, col]=size(H);
    out=cell(row,1);
    for i = 1:row
        out{i,1} = [];
    end
    for i=1:row
        for j=1:col
            if(H(i,j)==1)
                out{i,1} = [out{i,1} j];
            end
        end
    end
end

function out=vn_cn(H)
    [row, col]=size(H);
    out=cell(col,1);
    for i = 1:col
        out{i,1} = [];
    end
    for i=1:col
        for j=1:row
            if(H(j,i)==1)
                out{i,1} = [out{i,1} j];
            end
        end
    end
end

function cword = nrldpc_encode(B,z,msg)
    [m,n] = size(B);
    cword = zeros(1,n*z);
    cword(1:(n-m)*z) = msg;
    % Double-diagonal encoding
    temp = zeros(1,z);
    for i = 1:4 %row 1 to 4
        for j = 1:n-m %message columns
            temp = mod(temp + mul_sh(msg((j-1)*z+1:j*z),B(i,j)),2);
        end
    end
    if B(2,n-m+1) == -1
        p1_sh = B(3,n-m+1);
    else
        p1_sh = B(2,n-m+1);
    end
    cword((n-m)*z+1:(n-m+1)*z) = mul_sh(temp,z-p1_sh); %p1
    % Find p2, p3, p4
    for i = 1:3
        temp = zeros(1,z);
        for j = 1:n-m+i
            temp = mod(temp + mul_sh(cword((j-1)*z+1:j*z),B(i,j)),2);
        end
        cword((n-m+i)*z+1:(n-m+i+1)*z) = temp;
    end
    % Remaining parities
    for i = 5:m
        temp = zeros(1,z);
        for j = 1:n-m+4
            temp = mod(temp + mul_sh(cword((j-1)*z+1:j*z),B(i,j)),2);        
        end
        cword((n-m+i-1)*z+1:(n-m+i)*z) = temp;    
    end
end

function y = mul_sh(x,k)
    if(k==-1)
        y = zeros(1,length(x));
    else
        y = [x(k+1:end) x(1:k)];
    end
end