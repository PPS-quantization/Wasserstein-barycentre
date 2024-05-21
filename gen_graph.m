function A = gen_graph(n,type,d)

switch type
    case 'path' 
        M_path = [];
        M_path(1,:) = [1 0 zeros(1,n-2)];
        M_path(2,:) = [1 0 0 zeros(1,n-3)];
        M_path(n,:) = [zeros(1,n-2) 1 0 ];
        for i=2:n-2
            M_path(i+1,:) = circshift(M_path(i,:)',1)';
        end
        A = sparse(M_path);
        C = 1./max(sum(A)'*ones(1,size(A,1)),ones(size(A,1),1)*sum(A)).*A;
        C = C + diag(ones(size(A,1),1)-sum(C')');
        A = sparse(1/2*eye(size(A,1)) + 1/2*C);
    case 'bi_path' 
        M_path = [];
        M_path(1,:) = [0 1 zeros(1,n-2)];
        M_path(2,:) = [1 0 1 zeros(1,n-3)];
        M_path(n,:) = [zeros(1,n-2) 1 0 ];
        for i=2:n-2
            M_path(i+1,:) = circshift(M_path(i,:)',1)';
        end
        A = sparse(M_path);
        C = 1./max(sum(A)'*ones(1,size(A,1)),ones(size(A,1),1)*sum(A)).*A;
        C = C + diag(ones(size(A,1),1)-sum(C')');
        A = sparse(1/2*eye(size(A,1)) + 1/2*C);
    case 'grid2'
        M_path = [];
        M_path(1,:) = [0 1 zeros(1,n-2)];
        M_path(2,:) = [1 0 1 zeros(1,n-3)];
        M_path(n,:) = [zeros(1,n-2) 1 0 ];
        for i=2:n-2
            M_path(i+1,:) = circshift(M_path(i,:)',1)';
        end
        B = kron(M_path,eye(n))+kron(eye(n),M_path);
        C = 1./max(sum(B)'*ones(1,n^2),ones(n^2,1)*sum(B)).*B;
        C = C + diag(ones(n^2,1)-sum(C')');
        A = sparse(1/2*eye(n^2) + 1/2*C);
    case 'grid3'
        A1 = gen_graph(n,'grid2' );
        A2 = gen_graph(n,'bi_path' );
        B = kron(A1,eye(n))+kron(eye(n^2),A2);
        B(B>0) = 1;
        A = B./(sum(B)'*ones(1,n^3));
    case 'complete'
        A = 1/n*ones(n);
    case 'cycle'
        M_circ = [];
        M_circ(1,:) = [0 1 zeros(1,n-3) 1];
        for i=1:n-1
            M_circ(i+1,:) = circshift(M_circ(i,:)',1)';
        end
        A=sparse(1/2*eye(n) + 1/2*M_circ);
        A = sparse(A./(sum(A')'*ones(1,size(A,1))));
        B = adjacency(graph(A));
        C = 1./max(sum(B)'*ones(1,size(B,1)),ones(size(B,1),1)*sum(B)).*B;
        C = C + diag(ones(size(B,1),1)-sum(C')');
        A = sparse(1/2*eye(size(B,1)) + 1/2*C);
    case 'small-world'
        A = [];
        A(1,:) = [0 1 zeros(1,n-3) 1];
        for i=1:n-1
            A(i+1,:) = circshift(A(i,:)',1)';
        end
        ind = find(A==0);
        s_ind = size(ind,1);
        x = rand(s_ind,1);
        y = ones(s_ind,1);
        y(x>d) = 0;
        A(ind(y==1)) = 1;
        A = A+A';
        C = 1./max(sum(A)'*ones(1,size(A,1)),ones(size(A,1),1)*sum(A)).*A;
        C = C + diag(ones(size(A,1),1)-sum(C')');
        A = sparse(1/2*eye(size(A,1)) + 1/2*C);
    case 'torus-d'
        A1 = gen_graph(n,'cycle',0);
        for i=1:d-1
            A1 = kron(A1,eye(n))+kron(eye(size(A1,1)),gen_graph(n,'cycle',0));
        end
        A = sparse(A1./(sum(A1)'*ones(1,size(A1,1))));
    case 'complete_bin_tree'
        A1 = zeros(n*2+1);
        B = [zeros(n,1) kron(eye(n),[1 1])];
        A1(1:n,1:2*n+1) = B;
        A1 = [A1 + A1'];
        A1 = A1 + eye(2*n+1);
        A1 = A1./(sum(A1')'*ones(1,2*n+1));
        A =  sparse(A1);
    case 'dumbbell'
        A1 = ones(n);
        A2 = kron(eye(2),A1);
        A2(n,n+1)=1;
        A2(n+1,n)=1;
        C = 1./max(sum(A2)'*ones(1,size(A2,1)),ones(size(A2,1),1)*sum(A2)).*A2;
        C = C + diag(ones(size(A2,1),1)-sum(C')');
        A = sparse(1/2*eye(size(A2,1)) + 1/2*C);
    case 'star'
        % first node is in the center
        A = zeros(n);
        A(1,:) = [0 1/(n-1)*ones(1,n-1)];
        A(:,1) = [0 1/(n-1)*ones(1,n-1)]';
        A = sparse(A + diag(ones(n,1)-sum(A')'));
        A = sparse(1/2*eye(n) + 1/2*A);
    case 'star2'
        % first node is in the center
        A1 = gen_graph(n,'star');
        A = [A1 zeros(n); zeros(n) A1];
        A(1,n+1) = 1;
        A(n+1,1) = 1;
        A = A./(sum(A')'*ones(1,2*n));
    case 'bolas'
        A1 = gen_graph(n,'complete',0);
        A2 = gen_graph(n,'bi_path',0);
        A = [A1 zeros(n) zeros(n) ; zeros(n) A2 zeros(n); zeros(n) zeros(n) A1];
        A(n,n+1) = 1;
        A(n+1,n) = 1;
        A(2*n,2*n+1) = 1;
        A(2*n+1,2*n) = 1;
        A = sparse(A./(sum(A)'*ones(1,size(A,1))));
        
    case 'random_geometric'
        % extra parameter is the radius of minimum connection
        % It checks that the generated graphs is connected
        x_pos = rand(n,1);
        y_pos = rand(n,1);
        D = pdist2([x_pos y_pos],[x_pos y_pos]);
        A = ones(n);
        A(D>d) = 0;
        A = sparse(A./(sum(A)'*ones(1,size(A,1))));
        while abs(sum(eigs(A,2))-2)<1e-10
            x_pos = rand(n,1);
            y_pos = rand(n,1);
            D = pdist2([x_pos y_pos],[x_pos y_pos]);
            A = ones(n);
            A(D>d) = 0;
            A = sparse(A./(sum(A')'*ones(1,size(A,1))));
        end
    case 'erdos-renyi'
        A1 = spones(triu(sprand(n,n,d),1));
        A1 = A1 + A1';
        A1 = A1 + eye(n);
        % Now extracting the connectec component
        G = graph(A1);
        bins = conncomp(G);
        bincounts = histc(bins,1:max(bins));
        [val,ind] = max(bincounts);
        H = subgraph(G,find(bins==ind));
        B = adjacency(H);
        C = 1./max(sum(B)'*ones(1,size(B,1)),ones(size(B,1),1)*sum(B)).*B;
        C = C + diag(ones(size(B,1),1)-sum(C')');
        A = sparse(1/2*eye(size(B,1)) + 1/2*C);
        
%         % Tadpole graph (circle graph connected to path)
%         s = 5;    %number of nodes in the circle
%         M_tadp = zeros(n);
%         M_tadp(1:s-3,2:s) = [diag(1/2*ones(1,s-3)) zeros(s-3,1) [1/2;zeros(s-3-1,1)]];
%         M_tadp(s-2,s-1) = [1/3];
%         M_tadp(s-1,s:s+1) = [1/3 1/3];
%         M_tadp(s+1:n-1,s+2:n) = diag(1/2*ones(1,n-s-1));
%         M_tadp = M_tadp+M_tadp';
%         ss= 1-sum(M_tadp')';
%         M_tadp = M_tadp + diag(ss);
%         
%         lazy metropolis
%         N_tadp = 1/2*eye(n)+1/2*M_tadp;
%         
% % Bolas graph  (two connected graphs joined by path)
%         s = 3;    %number of nodes in the connected graphs
%         M_bolas = zeros(n);
%         M_bolas(1:s-2,2:s-1) = 1/(s-1)*triu(ones(s-2));
%         M_bolas(1:s-1,s) = 1/s*ones(s-1,1);
%         M_bolas(s:n-s,s+1:n-s+1) = diag([1/s 1/2*ones(1,n-2*s-1) 1/s]);
%         M_bolas(n-s+1,n-s+2:n)  = 1/s*ones(1,s-1);
%         M_bolas(n-s+2:n-1,n-s+3:n) = 1/(s-1)*triu(ones(s-2));
%         
%         M_bolas = M_bolas+M_bolas';
%         ss= 1-sum(M_bolas')';
%         M_bolas = M_bolas + diag(ss);
%         
%         N_bolas = 1/2*eye(n)+1/2*M_bolas;
% %% Lollipop graph   (connected graph joined by path)
%         s = 5;    %number of nodes in the connedted lollipop
%         M_loll = zeros(n);
%         M_loll(1,:) = [1/(s*(s-1)) 1/(s-1)*ones(1,s-2) 1/s zeros(1,n-s)];
%         for i=1:s-1
%             M_loll(i+1,:) = [circshift(M_loll(i,1:s-1)',1)' 1/s zeros(1,n-s)];
%         end
%         M_loll(s,:) = [1/s*ones(1,s-1) 0 1/s zeros(1,n-s-1)];
%         M_loll(s+1,:) = [zeros(1,s-1) 1/s (s-2)/(2*s) 1/2 zeros(1,n-s-2)];
%         M_loll(s+2,:) = [zeros(1,s) 1/2 0 1/2 zeros(1,n-s-3)];
%         for i=s+2:n-1
%             M_loll(i+1,:) = [zeros(1,s)  circshift(M_loll(i,s+1:n)',1)'];
%         end
%         M_loll(n,:) = [zeros(1,n-2) 1/2 1/2 ];
%         % lazy metropolis
%         N_loll = 1/2*eye(n)+1/2*M_loll;
    otherwise
        warning('Unexpected graph type. No graph created.')
end

