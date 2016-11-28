function Gadj = load_artif_graph(network_name)
%% artificial core network topology
% Ring Ring
if strcmpi('ringring', network_name)==true
    net = [
        1,2,1;
		1,5,1;
        2,3,1;
        3,4,1;
        3,7,1;
        4,5,1;
        4,30,1;
        6,7,1;
        6,10,1;
        7,8,1;
        7,21,1;
        8,9,1;
        8,14,1;
        9,10,1;
        11,12,1;
        11,14,1;
        12,13,1;
		12,16,1;
        12,21,1;
        13,15,1;
        14,15,1;
        16,17,1;
		16,20,1;
		17,18,1;
		18,19,1;
		19,20,1;
		19,24,1;
		19,28,1;
		21,22,1;
		21,25,1;
		22,23,1;
		23,24,1;
		24,25,1;
		24,30,1;
		26,27,1;
		26,30,1;
		27,28,1;
		28,29,1;
		29,30,1;
        ]+1;
else
    error('Invalid backbone name.');
end
num_node = max(max(net(:,1:2)));
Gadj = zeros(num_node);
for i = 1:size(net,1)
    Gadj(net(i,1),net(i,2)) = net(i,3);
end
Gadj = Gadj+Gadj';
end