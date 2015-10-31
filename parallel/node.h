#ifndef _NODE_
#define _NODE_

#include "boundary.h" // Any node has some boundaries with the neighboring nodes. Boundaries have information about adjasent nodes id and cells that are neighbor.

struct Node{
	int total_nodes; // total number of nodes
	int node_id; // node_id is the id of thisnode.
// Box is divided to reagions for our nodes. We lable the node position by idx and idy
	int idx,idy; // idx and idy shows the x position of the node inside the box
	int size_x, size_y; // size_x and size_y is the column and row number of cells within thisnode. It might be different, because the number of columns in the box is not dividable to the number of node columns.
// head_cell_idx is the idx of the first cell (the left most) in thisnode. The same for the head_cell_idy
// tial_cell_idx is the idx of the righ most cell plus one. The same for idy
	int head_cell_idx, head_cell_idy, tail_cell_idx, tail_cell_idy;
	long int seed; // seed number for initialization.
	int N; // Number of particles in the box. This will be transmitted from the box.
	Particle* particle; // This is a pointer to the original particle array pointer of the box. We need this pointer in some subroutins
	vector<Boundary> boundary; // Boundary list

// Static allocation (is not good for big simualations):
//Cell cell[divisor_x][divisor_y]; // We used cell list in our program. we divide the box to divisor_x by divisor_y cells. each cell has the information about particles id that are inside them.
// Dynamic allocation
	Cell** cell;
	
	Node();

	void Get_Box_Info(int size, Particle* p);
	void Init_Topology();
	void Send_Receive_Data(); // Send and Receive data of each neighboring cell
	void Quick_Update_Cells(); // Update particles that are inside each cell
	void Full_Update_Cells(); // Befor this function, Gather and Bcast must be called to have appropirate behaviour.
	void Update_Self_Neighbor_List(); // Updating neighborlist of particles inside cells within this node. But the pairs inside the node are considered
	void Update_Boundary_Neighbor_List(); // Updating neighborlist of particles inside cells within this node. But one the particles is outside this node.
	void Update_Neighbor_List(); // Updating neighborlist of particles inside cells within this node. All the pairs are considered.
	void Send_To_Root(); // Send thisnode information (particle position and angles) to the root node.
	void Root_Receive(); // Receive the sent information by other nodes
	void Root_Gather(); // Gather the information by root. Like a Send_To_Root() and Root_Receive() function.
	void Root_Bcast(); // Send all informations in root to other nodes

	void Neighbor_List_Interact(); // Interact using neighbor list
	void Self_Interact(); // Compute interaction of particles withing thisnode
	void Boundary_Interact(); // Compute interaction of particles of thisnode at the boundaries with the particles of neighboring node at the sam boundary.
	void Move(); // Move all particles within thisnode

	bool Chek_Seeds(); // Check if all nodes have seed number different from one another.
	void Print_Info(); // Print information of this node
};

Node::Node()
{
// Get the information abount total nodes and thisnode id
	MPI_Comm_size(MPI_COMM_WORLD, &total_nodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &node_id);

// Dynamical array allocation for cell:
	cell = new Cell*[divisor_x];
	for (int i = 0; i < divisor_x; i++)
		cell[i] = new Cell[divisor_y];

	for (int i = 0; i < divisor_x; i++)
		for (int j = 0; j < divisor_y; j++)
			cell[i][j].Init((Real) Lx*(2*i-divisor_x + 0.5)/divisor_x, (Real) Ly*(2*j-divisor_y + 0.5)/divisor_y); // setting the center position of each cell
}

void Node::Get_Box_Info(int size, Particle* p)
{
	N = size;
	particle = p;
	Cell::particle = p; // Each cell has a pointer to partilce array of the box. The cell needs this pointer for sum of its actions.
}

void Node::Init_Topology() // This function must be called after box definition.
{
// Check if the number of total nodes is in agreement with the way that the system is devided
	if ((npx*npy) != total_nodes)
	{
		cout << "Error, bad number of processors. Grid dimension is " << npx << " by " << npy << " that needs exactly " << npx*npy << " processors, but you ran the program with " << total_nodes << " processors." << endl;
		exit(0);
	}

// Computing the typical column and row number of cells in each node
	int width_x = divisor_x / npx;
	int remain_x = divisor_x % npx;
	int width_y = divisor_y / npy;
	int remain_y = divisor_y % npy;

// Finding the position of the node in the grid of nodes. First y index is increasing that means it changes faster than x index
	idx = node_id / npy;
	idy = node_id % npy;

// The first nodes are slightly bigger (width + 1) to match the size of the system. Their number is the same as the reminder of cells in division
	if (idx < remain_x)
	{
		head_cell_idx = idx*(width_x+1);
		size_x = width_x+1;
	}
	else
	{
		head_cell_idx = remain_x + idx*width_x;
		size_x = width_x;
	}

	if (idy < remain_y)
	{
		head_cell_idy = idy*(width_y+1);
		size_y = width_y+1;
	}
	else
	{
		head_cell_idy = remain_y + idy*width_y;
		size_y = width_y;
	}

// Finding the last cell idx and idy in the node
	tail_cell_idx = head_cell_idx + size_x;
	tail_cell_idy = head_cell_idy + size_y;

	int list_of_node[npx][npy]; // We need id of the other nodes by giving their position on grid
	for (int i = 0; i < npx; i++)
		for (int j = 0; j < npy; j++)
			list_of_node[i][j] = i*npy+j;


// Define a tmporary boundary object for pushback to boundary list of nodes
	Boundary temp_boundary;

// We first add all boundaries like a periodic boundary system, but at the end we will remove teh nodes boundaries that are within the box boundary.
// How we label boundaries are important. Typically (Like periodic boundary condition) each node must have 8 neighbors. Right, up righr, up, up left, left, down left, down, down right. We lable them as 0,1,2,3,4,5,6,7 respectively. Shecmaticly they will look like this:
//	3			2			1
//	4		thisnode		0
//	5			6			7


// 0 right edge
	temp_boundary.this_node_id = node_id; // node_id is the current node, thus we give it to the this node of the boundary.
	temp_boundary.that_node_id = list_of_node[(idx+1)%npx][idy]; // The right node id. Periodic transformation is applied.
	if (idx == npx)
		temp_boundary.box_edge = true;
	for (int y = head_cell_idy; y < tail_cell_idy; y++)
	{
		temp_boundary.this_cell.push_back(&cell[tail_cell_idx-1][y]); // Note: tail_cell_idx-1 is the id of the last cell in thisnode. Because tail_Cell_idx is greater than 1 we don't need % npx
		temp_boundary.that_cell.push_back(&cell[tail_cell_idx % divisor_x][y]); // Left node cell boundary x index is equal to tail_cell_idx
	}
	boundary.push_back(temp_boundary); // Push back to boundary
	temp_boundary.Delete(); // Deleting temporary boundary for next assignment
// 1 top right corner
	temp_boundary.this_node_id = node_id;
	temp_boundary.that_node_id = list_of_node[(idx+1)%npx][(idy+1)%npy];  // The above right node id. Periodic transformation is applied.
	if (idx == npx || idy == npy)
		temp_boundary.box_edge = true;
	temp_boundary.this_cell.push_back(&cell[tail_cell_idx-1][tail_cell_idy-1]);
	temp_boundary.that_cell.push_back(&cell[tail_cell_idx % divisor_x][tail_cell_idy % divisor_y]);
	boundary.push_back(temp_boundary);
	temp_boundary.Delete();
// 2 top edge
	temp_boundary.this_node_id = node_id;
	temp_boundary.that_node_id = list_of_node[idx][(idy+1)%npy]; // The above node id. Periodic transformation is applied.
	if (idy == npy)
		temp_boundary.box_edge = true;
	for (int x = head_cell_idx; x < tail_cell_idx; x++)
	{
		temp_boundary.this_cell.push_back(&cell[x][tail_cell_idy-1]);
		temp_boundary.that_cell.push_back(&cell[x][tail_cell_idy % divisor_y]);
	}
	boundary.push_back(temp_boundary);
	temp_boundary.Delete();
// 3 top left corner
	temp_boundary.this_node_id = node_id;
	temp_boundary.that_node_id = list_of_node[(idx-1+npx)%npx][(idy+1)%npy]; // The above left node id. Periodic transformation is applied.
	if (idx == 0 || idy == npy)
		temp_boundary.box_edge = true;
	temp_boundary.this_cell.push_back(&cell[head_cell_idx][tail_cell_idy-1]);
	temp_boundary.that_cell.push_back(&cell[(head_cell_idx-1+divisor_x) % divisor_x][tail_cell_idy % divisor_y]);
	boundary.push_back(temp_boundary);
	temp_boundary.Delete();
// 4 left edge
	temp_boundary.this_node_id = node_id;
	temp_boundary.that_node_id = list_of_node[(idx-1+npx)%npx][idy]; // The left node id. Periodic transformation is applied.
	if (idx == 0)
		temp_boundary.box_edge = true;
	for (int y = head_cell_idy; y < tail_cell_idy; y++)
	{
		temp_boundary.this_cell.push_back(&cell[head_cell_idx][y]);
		temp_boundary.that_cell.push_back(&cell[(head_cell_idx-1+divisor_x) % divisor_x][y]);
	}
	boundary.push_back(temp_boundary);
	temp_boundary.Delete();
// 5 bot left corner
	temp_boundary.this_node_id = node_id;
	temp_boundary.that_node_id = list_of_node[(idx-1+npx)%npx][(idy-1+npy)%npy]; // The down left node id. Periodic transformation is applied.
	if (idx == 0 || idy == 0)
		temp_boundary.box_edge = true;
	temp_boundary.this_cell.push_back(&cell[head_cell_idx][head_cell_idy]);
	temp_boundary.that_cell.push_back(&cell[(head_cell_idx-1+divisor_x) % divisor_x][(head_cell_idy-1+divisor_y) % divisor_y]);
	boundary.push_back(temp_boundary);
	temp_boundary.Delete();
// 6 bot edge
	temp_boundary.this_node_id = node_id;
	temp_boundary.that_node_id = list_of_node[idx][(idy-1+npy)%npy]; // The down node id. Periodic transformation is applied.
	if (idy == 0)
		temp_boundary.box_edge = true;
	for (int x = head_cell_idx; x < tail_cell_idx; x++)
	{
		temp_boundary.this_cell.push_back(&cell[x][head_cell_idy]);
		temp_boundary.that_cell.push_back(&cell[x][(head_cell_idy-1+divisor_y)%divisor_y]);
	}
	boundary.push_back(temp_boundary);
	temp_boundary.Delete();
// 7 bot right corner
	temp_boundary.this_node_id = node_id;
	temp_boundary.that_node_id = list_of_node[(idx+1)%npx][(idy-1+npy)%npy]; // The douwn right node id. Periodic transformation is applied.
	if (idx == npx || idy == 0)
		temp_boundary.box_edge = true;
	temp_boundary.this_cell.push_back(&cell[tail_cell_idx-1][head_cell_idy]);
	temp_boundary.that_cell.push_back(&cell[tail_cell_idx % divisor_x][(head_cell_idy-1+divisor_y) % divisor_y]);
	boundary.push_back(temp_boundary);
	temp_boundary.Delete();
// Assigning a uniqe tag to each boundary. The same boundary of different nodes must have the same tag.
	for (int i = 0; i < boundary.size(); i++)
	{
		if (boundary[i].this_node_id < boundary[i].that_node_id)
		{
			#ifdef PERIODIC_BOUNDARY_CONDITION
				boundary[i].tag = 0*boundary[i].box_edge*total_nodes*total_nodes + boundary[i].this_node_id*total_nodes + boundary[i].that_node_id;
			#else
				boundary[i].tag = boundary[i].this_node_id*total_nodes + boundary[i].that_node_id;
			#endif
		}
		else
		{
			#ifdef PERIODIC_BOUNDARY_CONDITION
				boundary[i].tag = 0*boundary[i].box_edge*total_nodes*total_nodes + boundary[i].that_node_id*total_nodes + boundary[i].this_node_id;
			#else
				boundary[i].tag = boundary[i].that_node_id*total_nodes + boundary[i].this_node_id;
			#endif
		}
	}
// Reforming boundaries in accord to the bondary condition.
	#ifndef PERIODIC_BOUNDARY_CONDITION
// Here the boundaries of edge nodes are removed:
	if (idx == (npx-1)) // The rightmost nodes
		if ((idy != 0) && (idy != npy-1)) // The nodes that are not at the corners (top and bot).
		{
// The erase must be from the last element, to avoid lable changes and a correct delletion for other boundaries
			boundary[7].is_active = false; // down right
			boundary[1].is_active = false; // up right
			boundary[0].is_active = false; // right
		}

	if (idx == 0) // The leftmost node
		if ((idy != 0) && (idy != npy-1)) // The nodes that are not at the corners (top and bot).
		{
			boundary[5].is_active = false; // down left
			boundary[4].is_active = false; // left
			boundary[3].is_active = false; // up left
		}

	if (idy == (npy-1)) // The top nodes
		if ((idx != 0) && (idx != npx-1)) // The nodes that are not at the corners (left and right).
		{
			boundary[3].is_active = false; // up left
			boundary[2].is_active = false; // up
			boundary[1].is_active = false; // up right
		}


	if (idy == 0) // The bottum nodes
		if ((idx != 0) && (idx != npx-1)) // The nodes that are not at the corners (left and right).
		{
			boundary[7].is_active = false; // down right
			boundary[6].is_active = false; // down
			boundary[5].is_active = false; // down left
		}

// Here the boundaries of corner nodes are removed:
	if (idx == 0) // The leftmost
		if (idy == 0) // The bottum
		{
			boundary[7].is_active = false;
			boundary[6].is_active = false;
			boundary[5].is_active = false;
			boundary[4].is_active = false;
			boundary[3].is_active = false;
		}
	if (idx == npx-1) // The rightmost
		if (idy == 0) // The bottum
		{
			boundary[7].is_active = false;
			boundary[6].is_active = false;
			boundary[5].is_active = false;
			boundary[1].is_active = false;
			boundary[0].is_active = false;
		}
	if (idx == npx-1) // The rightmost
		if (idy == npy-1) // The top
		{
			boundary[7].is_active = false;
			boundary[3].is_active = false;
			boundary[2].is_active = false;
			boundary[1].is_active = false;
			boundary[0].is_active = false;
		}
	if (idx == 0) // The leftmost
		if (idy == npy-1) // The top
		{
			boundary[5].is_active = false;
			boundary[4].is_active = false;
			boundary[3].is_active = false;
			boundary[2].is_active = false;
			boundary[1].is_active = false;
		}
	#endif
	MPI_Barrier(MPI_COMM_WORLD);
// All nodes are ready
}

// Send_Receive_Data will update boundary cells of each node with its neighboring nodes
// each node sends its information of boundary cells to the correspounding node.
void Node::Send_Receive_Data()
{
	for (int i = 0; i < boundary.size(); i++)
	{
		if (i % 4 != 2)
		{
			if (idx % 2 == 0)
			{
				if (boundary[i].is_active)
					boundary[i].Send_Data(); // Send information of i'th boundary of thisnode to the neighboring node that shares this boundary.
			}
			else
				if (boundary[(i+4)%8].is_active)
					boundary[(i+4)%8].Receive_Data(); // Receive information of i'th boundary of neighboring node.
			if (idx % 2 == 1)
			{
				if (boundary[i].is_active)
					boundary[i].Send_Data(); // Send information of i'th boundary of thisnode to the neighboring node that shares this boundary (if there is any).
			}
			else
				if (boundary[(i+4)%8].is_active)
					boundary[(i+4)%8].Receive_Data(); // Receive information of i'th boundary of neighboring node (if the neighbor exitst).
		}
		else
		{
			if (idy % 2 == 0)
			{
				if (boundary[i].is_active)
					boundary[i].Send_Data(); // Send information of i'th boundary of thisnode to the neighboring node that shares this boundary (if there is any).
			}
			else
				if (boundary[(i+4)%8].is_active)
					boundary[(i+4)%8].Receive_Data(); // Receive information of i'th boundary of neighboring node (if the neighbor exitst).
			if (idy % 2 == 1)
			{
				if (boundary[i].is_active)
					boundary[i].Send_Data(); // Send information of i'th boundary of thisnode to the neighboring node that shares this boundary (if there is any).
			}
			else
				if (boundary[(i+4)%8].is_active)
					boundary[(i+4)%8].Receive_Data(); // Receive information of i'th boundary of neighboring node (if the neighbor exitst).
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

// Quick_Update_Cells will update cells of each node (their particle) with the local information that means we have only information about particle position of thisnode and the boundary cells. This must be quicker than usage of the global information with a gather and bcast. Here neighobr list of each particle is computed as well.
void Node::Quick_Update_Cells()
{
	Send_Receive_Data();
	MPI_Barrier(MPI_COMM_WORLD);

	vector<int> node_pid; // pid is particle ids that possibly are within this node

// First we add particles in the neighboring cells which are not within the node. These particles may travell inside thisnode and we add them to the list of possible particles (node_pid). thisnode has a list of boundaries (right, top right, ...) and in the list of boundaries we have pointer to cells that belong to thisnode (this_cell) or to the neighobring node (that_cell). Here we only add that_cell particle ids because in future we will add all thisnode particles.
	for (int i = 0; i < boundary.size(); i++)
		for (int j = 0; j < boundary[i].that_cell.size(); j++)
		{
			for (int k = 0; k < boundary[i].that_cell[j]->pid.size(); k++)
				node_pid.push_back(boundary[i].that_cell[j]->pid[k]);
		}
	for (int i = 0; i < boundary.size(); i++)
		for (int j = 0; j < boundary[i].that_cell.size(); j++)
			boundary[i].that_cell[j]->Delete();
			

// In this part we go over all cells of thisnode (the first two for) and add particle of each cell to the node_pid (particle ids that may be inside thisnode). Here I also delete each cell, because I will add particles in the node_pid to the cell the cell must be empty to avoid repeatation of particles list in a cell.
	for (int x = head_cell_idx; x < tail_cell_idx; x++)
		for (int y = head_cell_idy; y < tail_cell_idy; y++)
		{
			for (int k = 0; k < cell[x][y].pid.size(); k++)
				node_pid.push_back(cell[x][y].pid[k]);
			cell[x][y].Delete();
		}

// Here the program checks each particle in node_pid. If the particle position is in a cell which belongs to thisnode, the program will add them to the list. It is very important that information about particles of neighboring cells must be up to date. For example this function must be used after an interaction computation to make sure that recently such an update has been occured.
	for (int i = 0; i < node_pid.size(); i++)
	{
// Find the index of the cell in which a particle are located.
		int x,y;
		x = (int) (particle[node_pid[i]].r.x + Lx)*divisor_x / Lx2;
		y = (int) (particle[node_pid[i]].r.y + Ly)*divisor_y / Ly2;

// Check if the particles are inside the box for a debug.
		#ifdef DEBUG
		if ((x >= divisor_x) || (x < 0) || (y >= divisor_y) || (y < 0))
		{
			cout << "\n Particle number " << node_pid[i] << " is Out of the box" << endl << flush;
			cout << "Particle Position is " << particle[node_pid[i]].r << endl;
			cout << "Particle  " << particle[node_pid[i]].v << endl;
			exit(0);
		}
		#endif

		#ifdef TRACK_PARTICLE
//			if (node_pid[i] == track)
//				if ((head_cell_idx <= x) && (x < tail_cell_idx))
//					if ((head_cell_idy <=  y)&& (y < tail_cell_idy))
//						if (flag)
//							cout << "Node: " << node_id << " Cell: " << x << " " << y << " " << particle[track].r << " " << particle[track].theta << endl << flush;
		#endif
		
		cell[x][y].Add(node_pid[i]);
	}

// We don't need node_pid anymore and we need it to be empty for our further use.
	node_pid.clear();


// Now particle indices are changed and we have to update information of boundaries. The particles of other nodes that are at boundaries
	MPI_Barrier(MPI_COMM_WORLD);

	for (int i = 0; i < boundary.size(); i++)
	{
		if (i % 4 != 2)
		{
			if (idx % 2 == 0)
			{
				if (boundary[i].is_active)
					boundary[i].Send_Particle_Ids(); // Send information of i'th boundary of thisnode to the neighboring node that shares this boundary.
			}
			else
				if (boundary[(i+4)%8].is_active)
					boundary[(i+4)%8].Receive_Particle_Ids(); // Receive information of i'th boundary of neighboring node.
			if (idx % 2 == 1)
			{
				if (boundary[i].is_active)
					boundary[i].Send_Particle_Ids(); // Send information of i'th boundary of thisnode to the neighboring node that shares this boundary (if there is any).
			}
			else
				if (boundary[(i+4)%8].is_active)
					boundary[(i+4)%8].Receive_Particle_Ids(); // Receive information of i'th boundary of neighboring node (if the neighbor exitst).
		}
		else
		{
			if (idy % 2 == 0)
			{
				if (boundary[i].is_active)
					boundary[i].Send_Particle_Ids(); // Send information of i'th boundary of thisnode to the neighboring node that shares this boundary (if there is any).
			}
			else
				if (boundary[(i+4)%8].is_active)
					boundary[(i+4)%8].Receive_Particle_Ids(); // Receive information of i'th boundary of neighboring node (if the neighbor exitst).
			if (idy % 2 == 1)
			{
				if (boundary[i].is_active)
					boundary[i].Send_Particle_Ids(); // Send information of i'th boundary of thisnode to the neighboring node that shares this boundary (if there is any).
			}
			else
				if (boundary[(i+4)%8].is_active)
					boundary[(i+4)%8].Receive_Particle_Ids(); // Receive information of i'th boundary of neighboring node (if the neighbor exitst).
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
}

// Full_Update_Cells will update cells of each node (their particle) with the global information that means the master node will gather information of all other nodes and broadcast the whole information to every nodes. Therefor each node has the information of any other node and is aware of all particles. After we check all particles to see to which cell they belong.
void Node::Full_Update_Cells()
{
// First we need to empty the cells from our particle ids. (To avoid degeneracies!)
	for (int x = 0; x < divisor_x; x++)
		for (int y = 0; y < divisor_y; y++)
			cell[x][y].Delete();

	for (int i = 0; i < N; i++)
	{
// Find the index of the cell in which a particle are located.
		int x,y;
		x = (int) (particle[i].r.x + Lx)*divisor_x / Lx2;
		y = (int) (particle[i].r.y + Ly)*divisor_y / Ly2;

// Check if the particles are inside the box for a debug.
		#ifdef DEBUG
		if ((x >= divisor_x) || (x < 0) || (y >= divisor_y) || (y < 0))
		{
			cout << "\n Particle number " << i << " is Out of the box" << endl << flush;
			cout << "Particle Position is " << particle[i].r << endl;
			exit(0);
		}
		#endif

		cell[x][y].Add(i);
	}
}

// Using the information of particles we update a list for each particle showing the neighboring particles. But we are considering the third newton law. That means particles within the same node are counted once as neighbor in the neighbor list of one of the two particles.
void Node::Update_Self_Neighbor_List()
{
// Each cell must interact with itself and 4 of its 8 neihbors that are right cell, up cell, righ up and right down. Because each intertion compute the torque to both particles we need to use 4 of the 8 directions.

// Self interaction
	for (int x = head_cell_idx; x < tail_cell_idx; x++)
		for (int y = head_cell_idy; y < tail_cell_idy; y++)
			cell[x][y].Neighbor_List();
// right, up and up right cells:
// The righmost cells and top cells must be excluded to avoid nieghbor node interactions.
	for (int x = head_cell_idx; x < (tail_cell_idx-1); x++)
		for (int y = head_cell_idy; y < (tail_cell_idy-1); y++)
		{
// No need for % divisor_x and % divisor_y because we are only considering interaction of cells inside thisnode.
			cell[x][y].Neighbor_List(&cell[x+1][y]);
			cell[x][y].Neighbor_List(&cell[x][y+1]);
			cell[x][y].Neighbor_List(&cell[x+1][y+1]);
		}

// The rest that I forgot:
	for (int x = head_cell_idx; x < (tail_cell_idx-1); x++)
		cell[x][tail_cell_idy-1].Neighbor_List(&cell[x+1][tail_cell_idy-1]);
	for (int y = head_cell_idy; y < (tail_cell_idy-1); y++)
		cell[tail_cell_idx-1][y].Neighbor_List(&cell[tail_cell_idx-1][y+1]);

// right down cell:
// The righmost cells and buttom cells must be excluded to avoid nieghbor node interactions.
	for (int x = head_cell_idx; x < (tail_cell_idx-1); x++)
		for (int y = head_cell_idy+1; y < (tail_cell_idy); y++)
			cell[x][y].Neighbor_List(&cell[x+1][y-1]);	
}

// Interaction of thisnode particles with particles outside of thisnode
void Node::Update_Boundary_Neighbor_List()
{
	#ifdef PERIODIC_BOUNDARY_CONDITION
// The first and last columns are excluded to avoid multiple interaction for the same pair of cells
	for (int x = (head_cell_idx+1); x < (tail_cell_idx-1); x++)
	{
// Buttom cells of thisnode interacting
		cell[x][head_cell_idy].Neighbor_List(&cell[x][(head_cell_idy-1+divisor_y)%divisor_y]);
		cell[x][head_cell_idy].Neighbor_List(&cell[(x+1)%divisor_x][(head_cell_idy-1+divisor_y)%divisor_y]);
		cell[x][head_cell_idy].Neighbor_List(&cell[(x-1+divisor_x)%divisor_x][(head_cell_idy-1+divisor_y)%divisor_y]);

// Top cells of thisnode interacting
		cell[x][tail_cell_idy-1].Neighbor_List(&cell[x][tail_cell_idy%divisor_y]);
		cell[x][tail_cell_idy-1].Neighbor_List(&cell[(x+1)%divisor_x][tail_cell_idy%divisor_y]);
		cell[x][tail_cell_idy-1].Neighbor_List(&cell[(x-1+divisor_x)%divisor_x][tail_cell_idy%divisor_y]);
	}

// The first and last rows are excluded to avoid multiple interaction for the same pair of cells
	for (int y = (head_cell_idy+1); y < (tail_cell_idy-1); y++)
	{
// Left cells of this node
		cell[head_cell_idx][y].Neighbor_List(&cell[(head_cell_idx-1+divisor_x) % divisor_x][y]);
		cell[head_cell_idx][y].Neighbor_List(&cell[(head_cell_idx-1+divisor_x) % divisor_x][(y+1)%divisor_y]);
		cell[head_cell_idx][y].Neighbor_List(&cell[(head_cell_idx-1+divisor_x) % divisor_x][(y-1+divisor_y)%divisor_y]);

// Right cells of this node
		cell[tail_cell_idx-1][y].Neighbor_List(&cell[tail_cell_idx % divisor_x][y]);
		cell[tail_cell_idx-1][y].Neighbor_List(&cell[tail_cell_idx % divisor_x][(y+1)%divisor_y]);
		cell[tail_cell_idx-1][y].Neighbor_List(&cell[tail_cell_idx % divisor_x][(y-1+divisor_y)%divisor_y]);
	}

// Interaction of corners:
// Left Buttom:
	cell[head_cell_idx][head_cell_idy].Neighbor_List(&cell[(head_cell_idx+1)%divisor_x][(head_cell_idy-1+divisor_y)%divisor_y]); // 7
	cell[head_cell_idx][head_cell_idy].Neighbor_List(&cell[head_cell_idx][(head_cell_idy-1+divisor_y)%divisor_y]); // 6
	cell[head_cell_idx][head_cell_idy].Neighbor_List(&cell[(head_cell_idx-1+divisor_x) % divisor_x][(head_cell_idy-1+divisor_y)%divisor_y]); // 5
	cell[head_cell_idx][head_cell_idy].Neighbor_List(&cell[(head_cell_idx-1+divisor_x) % divisor_x][head_cell_idy]); // 4
	cell[head_cell_idx][head_cell_idy].Neighbor_List(&cell[(head_cell_idx-1+divisor_x) % divisor_x][(head_cell_idy+1)%divisor_y]); // 3
// Left Top:
	cell[head_cell_idx][(tail_cell_idy-1+divisor_y)%divisor_y].Neighbor_List(&cell[(head_cell_idx-1+divisor_x) % divisor_x][(tail_cell_idy-2+divisor_y)%divisor_y]); // 5
	cell[head_cell_idx][(tail_cell_idy-1+divisor_y)%divisor_y].Neighbor_List(&cell[(head_cell_idx-1+divisor_x) % divisor_x][(tail_cell_idy-1+divisor_y)%divisor_y]); // 4
	cell[head_cell_idx][(tail_cell_idy-1+divisor_y)%divisor_y].Neighbor_List(&cell[(head_cell_idx-1+divisor_x) % divisor_x][tail_cell_idy%divisor_y]); // 3
	cell[head_cell_idx][(tail_cell_idy-1+divisor_y)%divisor_y].Neighbor_List(&cell[head_cell_idx][tail_cell_idy%divisor_y]); // 2
	cell[head_cell_idx][(tail_cell_idy-1+divisor_y)%divisor_y].Neighbor_List(&cell[(head_cell_idx+1)%divisor_x][tail_cell_idy%divisor_y]); // 1
// Right Top:
	cell[tail_cell_idx-1][tail_cell_idy-1].Neighbor_List(&cell[(tail_cell_idx-2+divisor_x)%divisor_x][tail_cell_idy%divisor_y]); // 3
	cell[tail_cell_idx-1][tail_cell_idy-1].Neighbor_List(&cell[tail_cell_idx-1][tail_cell_idy%divisor_y]); // 2
	cell[tail_cell_idx-1][tail_cell_idy-1].Neighbor_List(&cell[tail_cell_idx % divisor_x][tail_cell_idy%divisor_y]); // 1
	cell[tail_cell_idx-1][tail_cell_idy-1].Neighbor_List(&cell[tail_cell_idx % divisor_x][tail_cell_idy-1]); // 0
	cell[tail_cell_idx-1][tail_cell_idy-1].Neighbor_List(&cell[tail_cell_idx % divisor_x][(tail_cell_idy-2+divisor_y)%divisor_y]); // 7
// Right Buttom:
	cell[tail_cell_idx-1][head_cell_idy].Neighbor_List(&cell[tail_cell_idx % divisor_x][(head_cell_idy+1)%divisor_y]); // 1
	cell[tail_cell_idx-1][head_cell_idy].Neighbor_List(&cell[tail_cell_idx % divisor_x][head_cell_idy]); // 0
	cell[tail_cell_idx-1][head_cell_idy].Neighbor_List(&cell[tail_cell_idx % divisor_x][(head_cell_idy-1+divisor_y)%divisor_y]); // 7
	cell[tail_cell_idx-1][head_cell_idy].Neighbor_List(&cell[tail_cell_idx-1][(head_cell_idy-1+divisor_y)%divisor_y]); // 6
	cell[tail_cell_idx-1][head_cell_idy].Neighbor_List(&cell[(tail_cell_idx-2+divisor_x)%divisor_x][(head_cell_idy-1+divisor_y)%divisor_y]); // 5
	#else
// In future I'm going to remove the if conditions for a better performance


// Right cells of thisnode interacting
	if (tail_cell_idx != divisor_x)
		for (int y = (head_cell_idy+1); y < (tail_cell_idy-1); y++) // The first and last rows are excluded to avoid multiple interaction for the same pair of cells
		{
			cell[tail_cell_idx-1][y].Neighbor_List(&cell[tail_cell_idx][y]);
			cell[tail_cell_idx-1][y].Neighbor_List(&cell[tail_cell_idx][y+1]);
			cell[tail_cell_idx-1][y].Neighbor_List(&cell[tail_cell_idx][y-1]);
		}

// Top cells of thisnode interacting
	if (tail_cell_idy != divisor_y) // If tail_cell y cordinate is not at the top
		for (int x = (head_cell_idx+1); x < (tail_cell_idx-1); x++) // The first and last columns are excluded to avoid multiple interaction for the same pair of cells
		{
			cell[x][tail_cell_idy-1].Neighbor_List(&cell[x][tail_cell_idy]);
			cell[x][tail_cell_idy-1].Neighbor_List(&cell[x+1][tail_cell_idy]);
			cell[x][tail_cell_idy-1].Neighbor_List(&cell[x-1][tail_cell_idy]);
		}

// Left cells of thisnode interacting
	if (head_cell_idx != 0)
		for (int y = (head_cell_idy+1); y < (tail_cell_idy-1); y++) // The first and last rows are excluded to avoid multiple interaction for the same pair of cells
		{
			cell[head_cell_idx][y].Neighbor_List(&cell[head_cell_idx-1][y]);
			cell[head_cell_idx][y].Neighbor_List(&cell[head_cell_idx-1][y+1]);
			cell[head_cell_idx][y].Neighbor_List(&cell[head_cell_idx-1][y-1]);
		}

// Buttom cells of thisnode interacting
	if (head_cell_idy != 0) // If head_cell y cordinate is not at the bottom
		for (int x = (head_cell_idx+1); x < (tail_cell_idx-1); x++) // The first and last columns are excluded to avoid multiple interaction for the same pair of cells
		{
			cell[x][head_cell_idy].Neighbor_List(&cell[x][head_cell_idy-1]);
			cell[x][head_cell_idy].Neighbor_List(&cell[x+1][head_cell_idy-1]);
			cell[x][head_cell_idy].Neighbor_List(&cell[x-1][head_cell_idy-1]);
		}

// Interaction of corners:
// Left Buttom:
	if (head_cell_idx != 0)
	{
		cell[head_cell_idx][head_cell_idy].Neighbor_List(&cell[head_cell_idx-1][head_cell_idy]); // left: Number 4
		cell[head_cell_idx][head_cell_idy].Neighbor_List(&cell[head_cell_idx-1][head_cell_idy+1]); // up left: Number 3
		if (head_cell_idy != 0)
			cell[head_cell_idx][head_cell_idy].Neighbor_List(&cell[head_cell_idx-1][head_cell_idy-1]); // down left: Number 5
	}
	if (head_cell_idy != 0)
	{
		cell[head_cell_idx][head_cell_idy].Neighbor_List(&cell[head_cell_idx][head_cell_idy-1]); // down: Number 6
		cell[head_cell_idx][head_cell_idy].Neighbor_List(&cell[head_cell_idx+1][head_cell_idy-1]); // down right: Number 7
	}
	
// Left Top:
	if (head_cell_idx != 0)
	{
		cell[head_cell_idx][tail_cell_idy-1].Neighbor_List(&cell[head_cell_idx-1][tail_cell_idy-1]); // left: Number 4
		cell[head_cell_idx][tail_cell_idy-1].Neighbor_List(&cell[head_cell_idx-1][tail_cell_idy-2]); // down left: Number 5
		if (tail_cell_idy != divisor_y)
			cell[head_cell_idx][tail_cell_idy-1].Neighbor_List(&cell[head_cell_idx-1][tail_cell_idy]); // up left: Number 3
	}
	if (tail_cell_idy != divisor_y)
	{
		cell[head_cell_idx][tail_cell_idy-1].Neighbor_List(&cell[head_cell_idx][tail_cell_idy]); // up: Number 2
		cell[head_cell_idx][tail_cell_idy-1].Neighbor_List(&cell[head_cell_idx+1][tail_cell_idy]); // up right: Number 1
	}

// Right Top:
	if (tail_cell_idx != divisor_x)
	{
		cell[tail_cell_idx-1][tail_cell_idy-1].Neighbor_List(&cell[tail_cell_idx][tail_cell_idy-1]); // right: Number 0
		cell[tail_cell_idx-1][tail_cell_idy-1].Neighbor_List(&cell[tail_cell_idx][tail_cell_idy-2]); // down right: Number 7
		if (tail_cell_idy != divisor_y)
			cell[tail_cell_idx-1][tail_cell_idy-1].Neighbor_List(&cell[tail_cell_idx][tail_cell_idy]); // up right: Number 1
	}
	if (tail_cell_idy != divisor_y)
	{
		cell[tail_cell_idx-1][tail_cell_idy-1].Neighbor_List(&cell[tail_cell_idx-1][tail_cell_idy]); // up: Number 2
		cell[tail_cell_idx-1][tail_cell_idy-1].Neighbor_List(&cell[tail_cell_idx-2][tail_cell_idy]); // up left: Number 3
	}

// Right Buttom:
	if (tail_cell_idx != divisor_x)
	{
		cell[tail_cell_idx-1][head_cell_idy].Neighbor_List(&cell[tail_cell_idx][head_cell_idy]); // right: Number 0
		cell[tail_cell_idx-1][head_cell_idy].Neighbor_List(&cell[tail_cell_idx][head_cell_idy+1]); // up right: Number 1
		if (head_cell_idy != 0)
			cell[tail_cell_idx-1][head_cell_idy].Neighbor_List(&cell[tail_cell_idx][head_cell_idy-1]); // down right: Number 7
	}
	if (head_cell_idy != 0)
	{
		cell[tail_cell_idx-1][head_cell_idy].Neighbor_List(&cell[tail_cell_idx-1][head_cell_idy-1]); // down: Number 6
		cell[tail_cell_idx-1][head_cell_idy].Neighbor_List(&cell[tail_cell_idx-2][head_cell_idy-1]); // down: Number 5
	}
	#endif
}

// This function must be called after transfer of data between nodes.
void Node::Update_Neighbor_List()
{
	for (int x = head_cell_idx; x < tail_cell_idx; x++)
		for (int y = head_cell_idy; y < tail_cell_idy; y++)
			cell[x][y].Clear_Neighbor_List();
	Update_Self_Neighbor_List();
	Update_Boundary_Neighbor_List();
}


// Sending information to Master node
void Node::Send_To_Root()
{
	if (node_id != 0)
	{
// Finding the number of particles within thisnode
		int particle_count = 0;
		for (int x = head_cell_idx; x < tail_cell_idx; x++)
			for (int y = head_cell_idy; y < tail_cell_idy; y++)
				particle_count += cell[x][y].pid.size();

// Allocation
		int* index_buffer = new int[particle_count];
		double* data_buffer = new double[3*particle_count];

// Assigning the buffer arrays
		int counter = 0; // count particles. We can not use a shift very simply becasue we need to write both index_buffer and data_buffer.
		for (int x = head_cell_idx; x < tail_cell_idx; x++)
			for (int y = head_cell_idy; y < tail_cell_idy; y++)
			{
				for (int i = 0; i < cell[x][y].pid.size(); i++)
				{
					int index = cell[x][y].pid[i];
					index_buffer[counter] = index;
					data_buffer[3*counter] = particle[index].r.x;
					data_buffer[3*counter+1] = particle[index].r.y;
					data_buffer[3*counter+2] = particle[index].theta;
					counter++;
				}
			}
// tag_max is the maximum of the available tag value. tag_max-1 is for index and tag_max is for the data
		MPI_Send(index_buffer, particle_count, MPI_INT, 0, tag_max-1,MPI_COMM_WORLD);
		MPI_Send(data_buffer, 3*particle_count, MPI_DOUBLE, 0, tag_max,MPI_COMM_WORLD);

// Deallocation
		delete [] index_buffer;
		delete [] data_buffer;
	}
}

void Node::Root_Receive()
{
// What master node does:
	if (node_id == 0)
	{
		// define a buffer for receiving data
		int*	index_buffer = new int[N]; // Indices of particles. The master node have to know which particle information it is receiving to save them properly and correctly.
		double* data_buffer; // data of particles. position and angle.
		int count; // number of particle that each node sends to masternode.
		for (int i = 1; i < total_nodes; i++)
		{
			MPI_Status status;
			MPI_Recv(index_buffer,N,MPI_INT,i,tag_max-1,MPI_COMM_WORLD,&status); // receiving the indices.
			MPI_Get_count(&status, MPI_INT, &count); // Finding the number of indices that masternode received form node i.
			data_buffer = new double[3*count]; // Initialize array with length 3*counts (3 double for each particle)
			MPI_Recv(data_buffer,3*count,MPI_DOUBLE,i,tag_max,MPI_COMM_WORLD,&status); // receiving the data
// Update each particle in according to the data that is received.
			for (int j = 0; j < count; j++)
			{
// particle id is index_buffer[j] and we assign the data to that particle x,y and theta.
				particle[index_buffer[j]].r.x = data_buffer[3*j];
				particle[index_buffer[j]].r.y = data_buffer[3*j+1];
				particle[index_buffer[j]].theta = data_buffer[3*j+2];
				particle[index_buffer[j]].v.x = cos(particle[index_buffer[j]].theta); // Optimization required, computing every particle velocities is not a good idea. A first step is computing the velocity of particles that are within the node, not the one on the neighboring cells of that node.
				particle[index_buffer[j]].v.y = sin(particle[index_buffer[j]].theta); // Optimization required, computing every particle velocities is not a good idea. A first step is computing the velocty of particles that are within the node, not the one on the neighboring cells of that node.
			}
			delete [] data_buffer; // We don't need data_buffer and because in the next for step the count may change, we need to initilize another buffer with the proper size.
		}
		delete [] index_buffer; // deleting the index buffer.
	}
}

// With this fucntion master node will gather the information of particles of any other node. In processes like saving the trajectory this is requiered.
void Node::Root_Gather()
{
	// Any node (thisnode) except the master node, must send its information to root (master node).
	Send_To_Root(); // Sending information to master node.
	Root_Receive(); // Receiving information by master node

	MPI_Barrier(MPI_COMM_WORLD);
}

// Bcast send the information of every particles from the master node to other nodes. Perhaps befor a Bcast we may call Gather to have the correct information of all particles.
void Node::Root_Bcast()
{
	double* data_buffer = new double[3*N]; // Data buffer, three times of particle number N (x,y and theta)
// Master node collect partilces information into the data_buffer.
	if (node_id == 0)
	{
		for (int i = 0; i < N; i++)
		{
			data_buffer[3*i] = particle[i].r.x;
			data_buffer[3*i+1] = particle[i].r.y;
			data_buffer[3*i+2] = particle[i].theta;
		}
	}
// Broad casting to all nodes. The root node is 0.
	MPI_Bcast(data_buffer, 3*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
// Other nodes have to assign the received valuse to the particles. No index is needed because we sent the information of particles by their order.
	if (node_id != 0)
	{
		for (int i = 0; i < N; i++)
		{
			particle[i].r.x = data_buffer[3*i];
			particle[i].r.y = data_buffer[3*i+1];
			particle[i].theta = data_buffer[3*i+2];
			particle[i].v.x = cos(particle[i].theta); // Optimization required, computing every particle velocities is not a good idea. A first step is computing the velocity of particles that are within the node, not the one on the neighboring cells of that node.
			particle[i].v.y = sin(particle[i].theta); // Optimization required, computing every particle velocities is not a good idea. A first step is computing the velocity of particles that are within the node, not the one on the neighboring cells of that node.
		}
	}
	delete [] data_buffer;
	MPI_Barrier(MPI_COMM_WORLD); // We want to make sure that all the nodes have the same information at the end (finished their task).
}

// Interaction of all particles within thisnode
void Node::Neighbor_List_Interact()
{
// Self interaction
	for (int x = head_cell_idx; x < tail_cell_idx; x++)
		for (int y = head_cell_idy; y < tail_cell_idy; y++)
			cell[x][y].Interact();
}

// Interaction of all particles within thisnode
void Node::Self_Interact()
{
// Each cell must interact with itself and 4 of its 8 neihbors that are right cell, up cell, righ up and right down. Because each intertion compute the torque to both particles we need to use 4 of the 8 directions.

// Self interaction
	for (int x = head_cell_idx; x < tail_cell_idx; x++)
		for (int y = head_cell_idy; y < tail_cell_idy; y++)
			cell[x][y].Self_Interact();
// right, up and up right cells:
// The righmost cells and top cells must be excluded to avoid nieghbor node interactions.
	for (int x = head_cell_idx; x < (tail_cell_idx-1); x++)
		for (int y = head_cell_idy; y < (tail_cell_idy-1); y++)
		{
// No need for % divisor_x and % divisor_y because we are only considering interaction of cells inside thisnode.
			cell[x][y].Interact(&cell[x+1][y]);
			cell[x][y].Interact(&cell[x][y+1]);
			cell[x][y].Interact(&cell[x+1][y+1]);
		}

// The rest that I forgot:
	for (int x = head_cell_idx; x < (tail_cell_idx-1); x++)
		cell[x][tail_cell_idy-1].Interact(&cell[x+1][tail_cell_idy-1]);
	for (int y = head_cell_idy; y < (tail_cell_idy-1); y++)
		cell[tail_cell_idx-1][y].Interact(&cell[tail_cell_idx-1][y+1]);

// right down cell:
// The righmost cells and buttom cells must be excluded to avoid nieghbor node interactions.
	for (int x = head_cell_idx; x < (tail_cell_idx-1); x++)
		for (int y = head_cell_idy+1; y < (tail_cell_idy); y++)
			cell[x][y].Interact(&cell[x+1][y-1]);
}

// Interaction of thisnode particles with particles outside of thisnode
void Node::Boundary_Interact()
{
	#ifdef PERIODIC_BOUNDARY_CONDITION
// The first and last columns are excluded to avoid multiple interaction for the same pair of cells
	for (int x = (head_cell_idx+1); x < (tail_cell_idx-1); x++)
	{
// Buttom cells of thisnode interacting
		cell[x][head_cell_idy].Interact(&cell[x][(head_cell_idy-1+divisor_y)%divisor_y]);
		cell[x][head_cell_idy].Interact(&cell[(x+1)%divisor_x][(head_cell_idy-1+divisor_y)%divisor_y]);
		cell[x][head_cell_idy].Interact(&cell[(x-1+divisor_x)%divisor_x][(head_cell_idy-1+divisor_y)%divisor_y]);

// Top cells of thisnode interacting
		cell[x][tail_cell_idy-1].Interact(&cell[x][tail_cell_idy%divisor_y]);
		cell[x][tail_cell_idy-1].Interact(&cell[(x+1)%divisor_x][tail_cell_idy%divisor_y]);
		cell[x][tail_cell_idy-1].Interact(&cell[(x-1+divisor_x)%divisor_x][tail_cell_idy%divisor_y]);
	}

// The first and last rows are excluded to avoid multiple interaction for the same pair of cells
	for (int y = (head_cell_idy+1); y < (tail_cell_idy-1); y++)
	{
// Left cells of this node
		cell[head_cell_idx][y].Interact(&cell[(head_cell_idx-1+divisor_x) % divisor_x][y]);
		cell[head_cell_idx][y].Interact(&cell[(head_cell_idx-1+divisor_x) % divisor_x][(y+1)%divisor_y]);
		cell[head_cell_idx][y].Interact(&cell[(head_cell_idx-1+divisor_x) % divisor_x][(y-1+divisor_y)%divisor_y]);

// Right cells of this node
		cell[tail_cell_idx-1][y].Interact(&cell[tail_cell_idx % divisor_x][y]);
		cell[tail_cell_idx-1][y].Interact(&cell[tail_cell_idx % divisor_x][(y+1)%divisor_y]);
		cell[tail_cell_idx-1][y].Interact(&cell[tail_cell_idx % divisor_x][(y-1+divisor_y)%divisor_y]);
	}

// Interaction of corners:
// Left Buttom:
	cell[head_cell_idx][head_cell_idy].Interact(&cell[(head_cell_idx+1)%divisor_x][(head_cell_idy-1+divisor_y)%divisor_y]); // 7
	cell[head_cell_idx][head_cell_idy].Interact(&cell[head_cell_idx][(head_cell_idy-1+divisor_y)%divisor_y]); // 6
	cell[head_cell_idx][head_cell_idy].Interact(&cell[(head_cell_idx-1+divisor_x) % divisor_x][(head_cell_idy-1+divisor_y)%divisor_y]); // 5
	cell[head_cell_idx][head_cell_idy].Interact(&cell[(head_cell_idx-1+divisor_x) % divisor_x][head_cell_idy]); // 4
	cell[head_cell_idx][head_cell_idy].Interact(&cell[(head_cell_idx-1+divisor_x) % divisor_x][(head_cell_idy+1)%divisor_y]); // 3
// Left Top:
	cell[head_cell_idx][(tail_cell_idy-1+divisor_y)%divisor_y].Interact(&cell[(head_cell_idx-1+divisor_x) % divisor_x][(tail_cell_idy-2+divisor_y)%divisor_y]); // 5
	cell[head_cell_idx][(tail_cell_idy-1+divisor_y)%divisor_y].Interact(&cell[(head_cell_idx-1+divisor_x) % divisor_x][(tail_cell_idy-1+divisor_y)%divisor_y]); // 4
	cell[head_cell_idx][(tail_cell_idy-1+divisor_y)%divisor_y].Interact(&cell[(head_cell_idx-1+divisor_x) % divisor_x][tail_cell_idy%divisor_y]); // 3
	cell[head_cell_idx][(tail_cell_idy-1+divisor_y)%divisor_y].Interact(&cell[head_cell_idx][tail_cell_idy%divisor_y]); // 2
	cell[head_cell_idx][(tail_cell_idy-1+divisor_y)%divisor_y].Interact(&cell[(head_cell_idx+1)%divisor_x][tail_cell_idy%divisor_y]); // 1
// Right Top:
	cell[tail_cell_idx-1][tail_cell_idy-1].Interact(&cell[(tail_cell_idx-2+divisor_x)%divisor_x][tail_cell_idy%divisor_y]); // 3
	cell[tail_cell_idx-1][tail_cell_idy-1].Interact(&cell[tail_cell_idx-1][tail_cell_idy%divisor_y]); // 2
	cell[tail_cell_idx-1][tail_cell_idy-1].Interact(&cell[tail_cell_idx % divisor_x][tail_cell_idy%divisor_y]); // 1
	cell[tail_cell_idx-1][tail_cell_idy-1].Interact(&cell[tail_cell_idx % divisor_x][tail_cell_idy-1]); // 0
	cell[tail_cell_idx-1][tail_cell_idy-1].Interact(&cell[tail_cell_idx % divisor_x][(tail_cell_idy-2+divisor_y)%divisor_y]); // 7
// Right Buttom:
	cell[tail_cell_idx-1][head_cell_idy].Interact(&cell[tail_cell_idx % divisor_x][(head_cell_idy+1)%divisor_y]); // 1
	cell[tail_cell_idx-1][head_cell_idy].Interact(&cell[tail_cell_idx % divisor_x][head_cell_idy]); // 0
	cell[tail_cell_idx-1][head_cell_idy].Interact(&cell[tail_cell_idx % divisor_x][(head_cell_idy-1+divisor_y)%divisor_y]); // 7
	cell[tail_cell_idx-1][head_cell_idy].Interact(&cell[tail_cell_idx-1][(head_cell_idy-1+divisor_y)%divisor_y]); // 6
	cell[tail_cell_idx-1][head_cell_idy].Interact(&cell[(tail_cell_idx-2+divisor_x)%divisor_x][(head_cell_idy-1+divisor_y)%divisor_y]); // 5
	#else
// In future I'm going to remove the if conditions for a better performance


// Right cells of thisnode interacting
	if (tail_cell_idx != divisor_x)
		for (int y = (head_cell_idy+1); y < (tail_cell_idy-1); y++) // The first and last rows are excluded to avoid multiple interaction for the same pair of cells
		{
			cell[tail_cell_idx-1][y].Interact(&cell[tail_cell_idx][y]);
			cell[tail_cell_idx-1][y].Interact(&cell[tail_cell_idx][y+1]);
			cell[tail_cell_idx-1][y].Interact(&cell[tail_cell_idx][y-1]);
		}

// Top cells of thisnode interacting
	if (tail_cell_idy != divisor_y) // If tail_cell y cordinate is not at the top
		for (int x = (head_cell_idx+1); x < (tail_cell_idx-1); x++) // The first and last columns are excluded to avoid multiple interaction for the same pair of cells
		{
			cell[x][tail_cell_idy-1].Interact(&cell[x][tail_cell_idy]);
			cell[x][tail_cell_idy-1].Interact(&cell[x+1][tail_cell_idy]);
			cell[x][tail_cell_idy-1].Interact(&cell[x-1][tail_cell_idy]);
		}

// Left cells of thisnode interacting
	if (head_cell_idx != 0)
		for (int y = (head_cell_idy+1); y < (tail_cell_idy-1); y++) // The first and last rows are excluded to avoid multiple interaction for the same pair of cells
		{
			cell[head_cell_idx][y].Interact(&cell[head_cell_idx-1][y]);
			cell[head_cell_idx][y].Interact(&cell[head_cell_idx-1][y+1]);
			cell[head_cell_idx][y].Interact(&cell[head_cell_idx-1][y-1]);
		}

// Buttom cells of thisnode interacting
	if (head_cell_idy != 0) // If head_cell y cordinate is not at the bottom
		for (int x = (head_cell_idx+1); x < (tail_cell_idx-1); x++) // The first and last columns are excluded to avoid multiple interaction for the same pair of cells
		{
			cell[x][head_cell_idy].Interact(&cell[x][head_cell_idy-1]);
			cell[x][head_cell_idy].Interact(&cell[x+1][head_cell_idy-1]);
			cell[x][head_cell_idy].Interact(&cell[x-1][head_cell_idy-1]);
		}

// Interaction of corners:
// Left Buttom:
	if (head_cell_idx != 0)
	{
		cell[head_cell_idx][head_cell_idy].Interact(&cell[head_cell_idx-1][head_cell_idy]); // left: Number 4
		cell[head_cell_idx][head_cell_idy].Interact(&cell[head_cell_idx-1][head_cell_idy+1]); // up left: Number 3
		if (head_cell_idy != 0)
			cell[head_cell_idx][head_cell_idy].Interact(&cell[head_cell_idx-1][head_cell_idy-1]); // down left: Number 5
	}
	if (head_cell_idy != 0)
	{
		cell[head_cell_idx][head_cell_idy].Interact(&cell[head_cell_idx][head_cell_idy-1]); // down: Number 6
		cell[head_cell_idx][head_cell_idy].Interact(&cell[head_cell_idx+1][head_cell_idy-1]); // down right: Number 7
	}
	
// Left Top:
	if (head_cell_idx != 0)
	{
		cell[head_cell_idx][tail_cell_idy-1].Interact(&cell[head_cell_idx-1][tail_cell_idy-1]); // left: Number 4
		cell[head_cell_idx][tail_cell_idy-1].Interact(&cell[head_cell_idx-1][tail_cell_idy-2]); // down left: Number 5
		if (tail_cell_idy != divisor_y)
			cell[head_cell_idx][tail_cell_idy-1].Interact(&cell[head_cell_idx-1][tail_cell_idy]); // up left: Number 3
	}
	if (tail_cell_idy != divisor_y)
	{
		cell[head_cell_idx][tail_cell_idy-1].Interact(&cell[head_cell_idx][tail_cell_idy]); // up: Number 2
		cell[head_cell_idx][tail_cell_idy-1].Interact(&cell[head_cell_idx+1][tail_cell_idy]); // up right: Number 1
	}

// Right Top:
	if (tail_cell_idx != divisor_x)
	{
		cell[tail_cell_idx-1][tail_cell_idy-1].Interact(&cell[tail_cell_idx][tail_cell_idy-1]); // right: Number 0
		cell[tail_cell_idx-1][tail_cell_idy-1].Interact(&cell[tail_cell_idx][tail_cell_idy-2]); // down right: Number 7
		if (tail_cell_idy != divisor_y)
			cell[tail_cell_idx-1][tail_cell_idy-1].Interact(&cell[tail_cell_idx][tail_cell_idy]); // up right: Number 1
	}
	if (tail_cell_idy != divisor_y)
	{
		cell[tail_cell_idx-1][tail_cell_idy-1].Interact(&cell[tail_cell_idx-1][tail_cell_idy]); // up: Number 2
		cell[tail_cell_idx-1][tail_cell_idy-1].Interact(&cell[tail_cell_idx-2][tail_cell_idy]); // up left: Number 3
	}

// Right Buttom:
	if (tail_cell_idx != divisor_x)
	{
		cell[tail_cell_idx-1][head_cell_idy].Interact(&cell[tail_cell_idx][head_cell_idy]); // right: Number 0
		cell[tail_cell_idx-1][head_cell_idy].Interact(&cell[tail_cell_idx][head_cell_idy+1]); // up right: Number 1
		if (head_cell_idy != 0)
			cell[tail_cell_idx-1][head_cell_idy].Interact(&cell[tail_cell_idx][head_cell_idy-1]); // down right: Number 7
	}
	if (head_cell_idy != 0)
	{
		cell[tail_cell_idx-1][head_cell_idy].Interact(&cell[tail_cell_idx-1][head_cell_idy-1]); // down: Number 6
		cell[tail_cell_idx-1][head_cell_idy].Interact(&cell[tail_cell_idx-2][head_cell_idy-1]); // down: Number 5
	}
	#endif
}

// Moving particles within thisnode
void Node::Move()
{
	for (int x = head_cell_idx; x < tail_cell_idx; x++)
		for (int y = head_cell_idy; y < tail_cell_idy; y++)
				cell[x][y].Move();
}

bool Node::Chek_Seeds()
{
	long int s[total_nodes];
	bool b;
	int int_b;
	MPI_Status status;
	if (node_id == 0)
	{
		s[0] = seed;
		for (int i = 1; i < total_nodes; i++)
			MPI_Recv(&s[i],1,MPI_LONG_INT,i,1,MPI_COMM_WORLD, &status);
	}
	else
		MPI_Send(&seed,1,MPI_LONG_INT,0,1,MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if (node_id == 0)
	{
		b = true;
		for (int i = 0; i < total_nodes; i++)
			for (int j = i+1; j < total_nodes; j++)
				b = b && (s[i] != s[j]);
		int_b = b;
		for (int i = 1; i < total_nodes; i++)
			MPI_Send(&int_b,1,MPI_INT,i,1,MPI_COMM_WORLD);
	}
	else
	{
		MPI_Recv(&int_b,1,MPI_INT,0,1,MPI_COMM_WORLD, &status);
		if (int_b == 1)
			b = true;
		else
			b = false;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	return (b);
}


void Node::Print_Info()
{
	cout << "Node: " << node_id << " cell dimx from: " << head_cell_idx << " to " << tail_cell_idx - 1 << " dimy from: " << head_cell_idy << " to " << tail_cell_idy - 1 << " Num. of boundaries: " << boundary.size() << " boundary nodes are: " << boundary[2].that_node_id << endl << flush;

}

#endif
