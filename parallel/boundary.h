#ifndef _BOUNDARY_
#define _BOUNDARY_

#include "mpi.h"

struct Boundary{
// Any node has a list of boundaries. Each boundary is aware of the node that it belongs to (this_node_id) and the node that it is connecting this_node_id to (that_node_id).
	int this_node_id;
	int that_node_id;
	int tag; // Tag is used for the message tags of our MPI. We must give to each boundary a uniqe tag to avoid confilicts.
	bool is_active; // This gives information about the boundary state, whether it is an inactive boundary (no real data transformation) or an active boundary (data must be transferred). For example if the boundary condition is a bounded box, the boundaries at the edges of the box are inactive becasue there is no neighboring node beyond the boundary.
	bool box_edge; // This give information about the boundary that is at the edge of the box or not
	vector<Cell*> this_cell; // the cells at the boundary that are in the this_node
	vector<Cell*> that_cell; // the cells at the boundary that are in the that_node

	Boundary();
	Boundary(const Boundary& b); // Copy constructor, because we want to manipulate boundaries by a vector (pushback) we need a copy constructor.
	~Boundary(); // We need temporary boundary objects, because a boundary object has pointer we have to make sure that the allocated space is freed to avoid memmory leak.
	void Delete(); // Freeing memory that is used by a boundary object

	void Send_Data(); // Send the data of particles in this_cell of boundary of this_node to that_cell of boundary of that_node. We suppose that the particle indices of that_cells are known exactly and use this to enhance our computation.
	void Receive_Data(); // Receive data of particles inside that cell wich are inside this_cell of that_node. We suppose that the particle indices of that_cells are known exactly and use this to enhance our computation.
	void Send_Particle_Ids(); // Send particle ids that are inside this_cell to the neighboring node. Befor sending data each node must be aware of its neighboring node particles at boundary.
	void Receive_Particle_Ids(); // Receive particle ids that are inside that_cell of neighboring node to this cell. Befor sending data each node must be aware of its neighboring node particles at boundary.
	void Print_Info();
};

Boundary::Boundary()
{
	tag = -1;
	box_edge = false;
	is_active = true;
}

// Copy constructor
Boundary::Boundary(const Boundary& b)
{
	Delete();
	this_node_id = b.this_node_id;
	that_node_id = b.that_node_id;
	tag = b.tag;
	is_active = b.is_active;
	for (int i = 0; i < b.this_cell.size(); i++)
		this_cell.push_back(b.this_cell[i]);
	for (int i = 0; i < b.that_cell.size(); i++)
		that_cell.push_back(b.that_cell[i]);
}

Boundary::~Boundary()
{
	Delete();
}

void Boundary::Delete()
{
	this_cell.clear();
	that_cell.clear();
	tag = -1;
	this_node_id = -1;
	that_node_id = -1;
}

void Boundary::Send_Data()
{
// First we need to find the data size.
	int data_size = 0;
// Go over all bondary cells:
	for (int i = 0; i < this_cell.size(); i++)
		data_size += this_cell[i]->pid.size(); // Summing particles of each neighboring cell
	data_size *= dof; // each particle has 3 double values (x,y and theta).
	double* data_buffer = new double[data_size]; // Allocating buffer array

	int shift = 0; // We need to have a track of the last element of data_buffer that we wrote.
// Go over particle ids of each boundary cell
	for (int i = 0; i < this_cell.size(); i++)
	{
		for (int j = 0; j < this_cell[i]->pid.size(); j++)
		{
			int index = this_cell[i]->pid[j]; // This is just for convinience. Save pid of j'th particle in i'th this_cell to index.
// save the index'th particle data to data_buffer
			data_buffer[shift+dof*j] = this_cell[i]->particle[index].r.x; // The particles could be accessed through Cell class
			data_buffer[shift+dof*j+1] = this_cell[i]->particle[index].r.y;
			data_buffer[shift+dof*j+2] = this_cell[i]->particle[index].theta;
			#ifdef NonPeriodicCompute
				data_buffer[shift+dof*j+3] = this_cell[i]->particle[index].r_original.x;
				data_buffer[shift+dof*j+4] = this_cell[i]->particle[index].r_original.y;
			#endif
		}
		shift += dof*this_cell[i]->pid.size(); // the last element id must be added with amount of data that we added in the for loop.
	}
	MPI_Send(data_buffer,data_size,MPI_DOUBLE,that_node_id,tag,MPI_COMM_WORLD);
	delete [] data_buffer;
}

void Boundary::Receive_Data()
{
// First we need to find the data size that thisnode is receiving.
	int data_size = 0;
// We go over the cells of neighboring node at the boundary to sum the number of particles.
	for (int i = 0; i < that_cell.size(); i++)
		data_size += that_cell[i]->pid.size();
	data_size *= dof; // each particle has 3 double values (x,y and theta).
	double* data_buffer = new double[data_size]; // Allocating space
	MPI_Status status; // status is required in MPI_Recv call
	MPI_Recv(data_buffer,data_size,MPI_DOUBLE,that_node_id,tag,MPI_COMM_WORLD,&status); // receiving data
	int shift = 0; // We need to have a track of the last element of data_buffer that we wrote.
	for (int i = 0; i < that_cell.size(); i++)
	{
		for (int j = 0; j < that_cell[i]->pid.size(); j++)
		{
			int index = that_cell[i]->pid[j];
			that_cell[i]->particle[index].r.x = data_buffer[shift+dof*j];
			that_cell[i]->particle[index].r.y = data_buffer[shift+dof*j+1];
			that_cell[i]->particle[index].theta = data_buffer[shift+dof*j+2];
			#ifdef NonPeriodicCompute
				that_cell[i]->particle[index].r_original.x = data_buffer[shift+dof*j+3];
				that_cell[i]->particle[index].r_original.y = data_buffer[shift+dof*j+4];
			#endif

			that_cell[i]->particle[index].v.x = cos(that_cell[i]->particle[index].theta); // Optimization required, computing every particle velocities is not a good idea. A first step is computing the velocity of particles that are within the node, not the one on the neighboring cells of that node.
			that_cell[i]->particle[index].v.y = sin(that_cell[i]->particle[index].theta); // Optimization required, computing every particle velocities is not a good idea. A first step is computing the velocity of particles that are within the node, not the one on the neighboring cells of that node.
			that_cell[i]->particle[index].Reset(); // eperimental for debug, it seems that this is needed!
		}
		shift += dof*that_cell[i]->pid.size();
	}
	delete [] data_buffer;
}

void Boundary::Send_Particle_Ids()
{
// First we need to find the data size.
	int index_size = 0;
// Go over all bondary cells:
	for (int i = 0; i < this_cell.size(); i++)
		index_size += this_cell[i]->pid.size(); // Summing particles of each neighboring cell
	int* index_buffer = new int[index_size]; // Allocating buffer array
	int* cell_size = new int[this_cell.size()]; // Allocating a buffer of each cell particle count.

	int shift = 0; // We need to have a track of the last element of data_buffer that we wrote.
// Go over particle ids of each boundary cell
	for (int i = 0; i < this_cell.size(); i++)
	{
		for (int j = 0; j < this_cell[i]->pid.size(); j++)
			index_buffer[shift+j] = this_cell[i]->pid[j];
		cell_size[i] = this_cell[i]->pid.size(); // Saving particle number of i'th cell of thisnode at boundary.
		shift += this_cell[i]->pid.size(); // the last element id must be added with amount of data that we added in the for loop.
	}
	MPI_Send(cell_size,this_cell.size(),MPI_INT,that_node_id,tag,MPI_COMM_WORLD);
	MPI_Send(index_buffer,index_size,MPI_INT,that_node_id,tag,MPI_COMM_WORLD);
	delete [] index_buffer;
	delete [] cell_size;
}

void Boundary::Receive_Particle_Ids()
{
	int* cell_size = new int[that_cell.size()]; // Allocating a buffer of each cell particle count.

	MPI_Status status; // status is required in MPI_Recv call
	MPI_Recv(cell_size, that_cell.size(),MPI_INT,that_node_id,tag,MPI_COMM_WORLD,&status); // Receiving particle number within each boundary cell

	int data_size = 0;
	for (int i = 0; i < that_cell.size(); i++)
		data_size += cell_size[i];
	int* index_buffer = new int[data_size]; // Allocating space
	MPI_Recv(index_buffer,data_size,MPI_DOUBLE,that_node_id,tag,MPI_COMM_WORLD,&status); // Receiving Indices

	int shift = 0; // We need to have a track of the last element of data_buffer that we wrote.
	for (int i = 0; i < that_cell.size(); i++)
	{
		that_cell[i]->Delete(); // Delete old informatinos
		for (int j = 0; j < cell_size[i]; j++)
		{
			int index = index_buffer[shift+j];
			that_cell[i]->pid.push_back(index);
		}
		shift += that_cell[i]->pid.size();
	}

	delete [] cell_size;
	delete [] index_buffer;
}

void Boundary::Print_Info()
{
	for (int i = 0; i < that_cell.size(); i++)
		for (int j = 0; j < that_cell[i]->pid.size(); j++)
			cout << "Boundary between node: " << this_node_id << " and " << that_node_id << " Particles are: " << that_cell[i]->pid[j] << endl << flush;
}

#endif
