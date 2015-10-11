/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012 Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include <cassert>
#include <algorithm>

#include "npdm_permutations.h"
#include "spinblock.h"
#include "pario.h"

namespace SpinAdapted{

//===========================================================================================================================================================

void Npdm_permutations::get_permute(const std::pair<std::vector<int>,int>& origin, int start, int n, std::vector<std::pair<std::vector<int>,int> >& reorders)
{
  if (n<2 ||start+1 >= n){
    reorders.push_back(origin);
    return;
  }
  get_permute(origin,start+1,n,reorders);
  for(int i=start+1; i<n;i++)
  {
    std::pair<std::vector<int>,int> neworder = origin;
    int tmp = neworder.first[i];
    neworder.first[i]= neworder.first[start];
    neworder.first[start] = tmp;
    neworder.second *= -1;
    get_permute(neworder,start+1,n,reorders);
  }


}


void Onepdm_permutations::get_spatial_batch( const std::vector< std::pair< std::vector<int>, double > >& in, 
                                              std::vector< std::pair< std::vector<int>, double > >& spatial_perms )
{
  spatial_perms.clear();
  std::vector<int> spatial(2);

  if (in.size()==0)
    return ;
  spatial[0] = in[0].first[0]/2;
  spatial[1] = in[0].first[1]/2;
  double value = in[0].second + in[1].second;
  spatial_perms.push_back(std::make_pair(spatial,value));
  if (spatial[0] !=spatial[1] && dmrginp.doimplicitTranspose())
    {
      int tmp = spatial[0];
      spatial[0] = spatial[1];
      spatial[1] = tmp;
      spatial_perms.push_back(std::make_pair(spatial,value));
    }
}

//===========================================================================================================================================================

void Twopdm_permutations::get_spatial_batch( const std::vector< std::pair< std::vector<int>, double > >& in, 
                                              std::vector< std::pair< std::vector<int>, double > >& spatial_perms )
{
  spatial_perms.clear();
  std::vector<int> indices;
  std::vector< std::vector<int> > newindices;
  std::vector<int> tmp;

  if (in.size()==0)
    return ;

  for(int i=0;i<in[0].first.size();i++)
    indices.push_back(in[0].first.at(i)/2);

  for(int i=0;i<reorders.size();i++)
  {
    
    std::vector<int> spatial_indices = indices;
    for(int j=0; j< reorders[i].first.size();j++)
    {
      spatial_indices[j] = indices[reorders[i].first[j]];
    }
    bool skip = false;
    for(int j=0; j< spatial_perms.size();j++)
    {
      if (spatial_perms[j].first == spatial_indices)
      {
        skip = true;
        break;
      }
    }
    if(!skip)
    {
      double value = 0;
      for(int j=0; j<in.size();j++)
      {
        if (in[j].first.at(reorders[i].first[0])%2 == in[j].first.at(3)%2 &&
            in[j].first.at(reorders[i].first[1])%2 == in[j].first.at(2)%2 )
          value += in[j].second*reorders[i].second;
      }
      if(abs(value)< NUMERICAL_ZERO)
        continue;
      std::set<std::vector<int> > spatial_batch;
      std::vector<int> tmp(4);
      for(int k=0;k<reorders.size();k++)
      {
        tmp[0]= spatial_indices[reorders[k].first[0]];
        tmp[1]= spatial_indices[reorders[k].first[1]];
        tmp[2]= spatial_indices[3-reorders[k].first[1]];
        tmp[3]= spatial_indices[3-reorders[k].first[0]];
        spatial_batch.insert(tmp);
        if(dmrginp.doimplicitTranspose())
        {
          tmp[0]= spatial_indices[3-reorders[k].first[0]];
          tmp[1]= spatial_indices[3-reorders[k].first[1]];
          tmp[2]= spatial_indices[reorders[k].first[1]];
          tmp[3]= spatial_indices[reorders[k].first[0]];
          spatial_batch.insert(tmp);
        }
      }
      for(std::set< std::vector<int> >::iterator x = spatial_batch.begin(); x != spatial_batch.end(); ++x)
        spatial_perms.push_back(std::make_pair(*x,value));
    }
  }
}

//===========================================================================================================================================================

void Threepdm_permutations::get_spatial_batch( const std::vector< std::pair< std::vector<int>, double > >& in, 
                                              std::vector< std::pair< std::vector<int>, double > >& spatial_perms )
{
  spatial_perms.clear();
  std::vector<int> indices;
  std::vector< std::vector<int> > newindices;
  std::vector<int> tmp;

  if (in.size()==0)
    return ;

  for(int i=0;i<in[0].first.size();i++)
    indices.push_back(in[0].first.at(i)/2);

  for(int i=0;i<reorders.size();i++)
  {
    
    std::vector<int> spatial_indices = indices;
    for(int j=0; j< reorders[i].first.size();j++)
    {
      spatial_indices[j] = indices[reorders[i].first[j]];
    }
    bool skip = false;
    for(int j=0; j< spatial_perms.size();j++)
    {
      if (spatial_perms[j].first == spatial_indices)
      {
        skip = true;
        break;
      }
    }
    if(!skip)
    {
      double value = 0;
      for(int j=0; j<in.size();j++)
      {
        if (in[j].first.at(reorders[i].first[0])%2 == in[j].first.at(5)%2 &&
            in[j].first.at(reorders[i].first[1])%2 == in[j].first.at(4)%2 &&
            in[j].first.at(reorders[i].first[2])%2 == in[j].first.at(3)%2 )
          value += in[j].second*reorders[i].second;
      }
      if(abs(value)< NUMERICAL_ZERO)
        continue;
      std::set<std::vector<int> > spatial_batch;
      std::vector<int> tmp(6);
      for(int k=0;k<reorders.size();k++)
      {
        tmp[0]= spatial_indices[reorders[k].first[0]];
        tmp[1]= spatial_indices[reorders[k].first[1]];
        tmp[2]= spatial_indices[reorders[k].first[2]];
        tmp[3]= spatial_indices[5-reorders[k].first[2]];
        tmp[4]= spatial_indices[5-reorders[k].first[1]];
        tmp[5]= spatial_indices[5-reorders[k].first[0]];
        spatial_batch.insert(tmp);
        if(dmrginp.doimplicitTranspose())
        {
          tmp[0]= spatial_indices[5-reorders[k].first[0]];
          tmp[1]= spatial_indices[5-reorders[k].first[1]];
          tmp[2]= spatial_indices[5-reorders[k].first[2]];
          tmp[3]= spatial_indices[reorders[k].first[2]];
          tmp[4]= spatial_indices[reorders[k].first[1]];
          tmp[5]= spatial_indices[reorders[k].first[0]];
          spatial_batch.insert(tmp);
        }
      }
   //   int k = spatial_indices[0];
   //   int l = spatial_indices[1];
   //   int m = spatial_indices[2];
   //   int n = spatial_indices[3];
   //   int p = spatial_indices[4];
   //   int q = spatial_indices[5];
   //   tmp={k,l,m,n,p,q}; spatial_batch.insert(tmp);
   //   tmp={k,m,l,p,n,q}; spatial_batch.insert(tmp);
   //   tmp={l,k,m,n,q,p}; spatial_batch.insert(tmp);
   //   tmp={l,m,k,q,n,p}; spatial_batch.insert(tmp);
   //   tmp={m,k,l,p,q,n}; spatial_batch.insert(tmp);
   //   tmp={m,l,k,q,p,n}; spatial_batch.insert(tmp);
   //   if(dmrginp.doimplicitTranspose())
   //   {
   //     tmp={q,p,n,m,l,k}; spatial_batch.insert(tmp);
   //     tmp={q,n,p,l,m,k}; spatial_batch.insert(tmp);
   //     tmp={p,q,n,m,k,l}; spatial_batch.insert(tmp);
   //     tmp={p,n,q,k,m,l}; spatial_batch.insert(tmp);
   //     tmp={n,q,p,l,k,m}; spatial_batch.insert(tmp);
   //     tmp={n,p,q,k,l,m}; spatial_batch.insert(tmp);
   //   }
      for(std::set< std::vector<int> >::iterator x = spatial_batch.begin(); x != spatial_batch.end(); ++x)
        spatial_perms.push_back(std::make_pair(*x,value));
    }
  }
}

//===========================================================================================================================================================

void Fourpdm_permutations::get_spatial_batch( const std::vector< std::pair< std::vector<int>, double > >& in, 
                                              std::vector< std::pair< std::vector<int>, double > >& spatial_perms )
{
  spatial_perms.clear();
  std::vector<int> indices;
  std::vector< std::vector<int> > newindices;
  std::vector<int> tmp;

  if (in.size()==0)
    return ;

  for(int i=0;i<in[0].first.size();i++)
    indices.push_back(in[0].first.at(i)/2);

  for(int i=0;i<reorders.size();i++)
  {
    
    std::vector<int> spatial_indices = indices;
    for(int j=0; j< reorders[i].first.size();j++)
    {
      spatial_indices[j] = indices[reorders[i].first[j]];
    }
    bool skip = false;
    for(int j=0; j< spatial_perms.size();j++)
    {
      if (spatial_perms[j].first == spatial_indices)
      {
        skip = true;
        break;
      }
    }
    if(!skip)
    {
      double value = 0;
      for(int j=0; j<in.size();j++)
      {
        if (in[j].first.at(reorders[i].first[0])%2 == in[j].first.at(7)%2 &&
            in[j].first.at(reorders[i].first[1])%2 == in[j].first.at(6)%2 &&
            in[j].first.at(reorders[i].first[2])%2 == in[j].first.at(5)%2 &&
            in[j].first.at(reorders[i].first[3])%2 == in[j].first.at(4)%2 )
          value += in[j].second*reorders[i].second;
      }
      if(abs(value)< NUMERICAL_ZERO)
        continue;
      std::vector<int> tmp(8);
      std::set<std::vector<int> > spatial_batch;
      for(int k=0;k<reorders.size();k++)
      {
        tmp[0]= spatial_indices[reorders[k].first[0]];
        tmp[1]= spatial_indices[reorders[k].first[1]];
        tmp[2]= spatial_indices[reorders[k].first[2]];
        tmp[3]= spatial_indices[reorders[k].first[3]];
        tmp[4]= spatial_indices[7-reorders[k].first[3]];
        tmp[5]= spatial_indices[7-reorders[k].first[2]];
        tmp[6]= spatial_indices[7-reorders[k].first[1]];
        tmp[7]= spatial_indices[7-reorders[k].first[0]];
        spatial_batch.insert(tmp);
        if(dmrginp.doimplicitTranspose())
        {
          tmp[0]= spatial_indices[7-reorders[k].first[0]];
          tmp[1]= spatial_indices[7-reorders[k].first[1]];
          tmp[2]= spatial_indices[7-reorders[k].first[2]];
          tmp[3]= spatial_indices[7-reorders[k].first[3]];
          tmp[4]= spatial_indices[reorders[k].first[3]];
          tmp[5]= spatial_indices[reorders[k].first[2]];
          tmp[6]= spatial_indices[reorders[k].first[1]];
          tmp[7]= spatial_indices[reorders[k].first[0]];
          spatial_batch.insert(tmp);
        }
      }
      for(std::set< std::vector<int> >::iterator x = spatial_batch.begin(); x != spatial_batch.end(); ++x)
        spatial_perms.push_back(std::make_pair(*x,value));
    }
  }
}

//===========================================================================================================================================================

void Npdm_permutations::process_new_elements( const std::vector< std::pair< std::vector<int>, double > >& in, 
                                              std::vector< std::pair< std::vector<int>, double > >& nonredundant_elements,
                                              std::vector< std::pair< std::vector<int>, double > >& spin_perms )
{
  spin_perms.clear();
  int count = 0;
   
  std::vector<int> spatial_indices;
  if (in.size() != 0)
  {
    for(int i=0;i<in[0].first.size();i++)
      spatial_indices.push_back(in[0].first.at(i)/2);

  }
  // Loop over all input spin indices
  for (int i=0; i<in.size(); i++) {
    // Get all permutations of each set of spin indices
    std::vector< std::pair< std::vector<int>, double > > tmp;
    get_spin_permutations( tmp, in[i].first, in[i].second );
    // If permutations do not generate any of the following input indices, save the non-redundant original
    bool keep = true;
    for ( int j=0; j<tmp.size(); j++) {
      // Loop over remaining original indices
      for (int k=i+1; k<in.size(); k++) {
        if ( tmp[j].first == in[k].first ) {
          keep = false;
          break;
        }
      }
      if (!keep) break;
    }
    if (keep) {
       count++;
       nonredundant_elements.push_back( in[i] );
       spin_perms.insert( spin_perms.end(), tmp.begin(), tmp.end() );
    }
  }

  assert( count <= in.size() );
//pout << "nonredundant elements = " << count << endl;

}

//===========================================================================================================================================================

void Onepdm_permutations::get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
                                                 const std::vector<int>& indices, const double& val )
{
  assert( indices.size() == 2 );
  std::vector<int> idx(2);
  int i = indices[0];
  int j = indices[1];

  idx[0] = i; idx[1] = j;
  spin_batch.push_back( std::make_pair( idx, val ) );
  // Transpose is same 
  if(dmrginp.doimplicitTranspose())
  if ( i != j ) {
    idx[0] = j; idx[1] = i;
    spin_batch.push_back( std::make_pair( idx, val ) );
  }

}

//===========================================================================================================================================================

void Twopdm_permutations::get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
                                                 const std::vector<int>& indices, const double& val )
{
  assert( indices.size() == 4 );
  std::vector<int> idx(4);
  int i = indices[0];
  int j = indices[1];
  int k = indices[2];
  int l = indices[3];

  // If indices are not all unique, then all elements should be zero (and next_even_permutation fails)
  std::vector<int> v(2);
  v[0] = i; v[1] = j;
  std::sort( v.begin(), v.end() );
  if ( (v[0]==v[1]) ) return;
  std::vector<int> w(2);
  w[0] = k; w[1] = l;
  std::sort( w.begin(), w.end() );
  if ( (w[0]==w[1]) ) return;
  bool skip_transpose = ( v == w );

  // 8 permutations
  //--------------------------
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; spin_batch.push_back( std::make_pair( idx, val ) );
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; spin_batch.push_back( std::make_pair( idx, val ) );
                    
  if ( !skip_transpose && dmrginp.doimplicitTranspose()) {
    idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; spin_batch.push_back( std::make_pair( idx, val ) );
    idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; spin_batch.push_back( std::make_pair( idx, val ) );
  }

}

//===========================================================================================================================================================

void Threepdm_permutations::get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
                                                   const std::vector<int>& indices, const double& val )
{
  assert( indices.size() == 6 );
  std::vector<int> idx(6);
  int i = indices[0];
  int j = indices[1];
  int k = indices[2];
  int l = indices[3];
  int m = indices[4];
  int n = indices[5];

  // If indices are not all unique, then all elements should be zero (and next_even_permutation fails)
  std::vector<int> v(3);
  v[0] = i; v[1] = j; v[2] = k;
  std::sort( v.begin(), v.end() );
  if ( (v[0]==v[1]) || (v[1]==v[2]) ) return;
  std::vector<int> w(3);
  w[0] = l; w[1] = m; w[2] = n;
  std::sort( w.begin(), w.end() );
  if ( (w[0]==w[1]) || (w[1]==w[2]) ) return;
  bool skip_transpose = ( v == w );

  // The number of possible permutations is (3!)**2 twice
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m; idx[5] = n; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = n; idx[5] = m; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = n; idx[4] = m; idx[5] = l; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = n; idx[4] = l; idx[5] = m; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = m; idx[4] = l; idx[5] = n; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = m; idx[4] = n; idx[5] = l; spin_batch.push_back( std::make_pair( idx,  val ) );

  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = m; idx[5] = n; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = n; idx[5] = m; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = n; idx[4] = m; idx[5] = l; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = n; idx[4] = l; idx[5] = m; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = m; idx[4] = l; idx[5] = n; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = m; idx[4] = n; idx[5] = l; spin_batch.push_back( std::make_pair( idx, -val ) );

  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = m; idx[5] = n; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = n; idx[5] = m; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = n; idx[4] = m; idx[5] = l; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = n; idx[4] = l; idx[5] = m; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = m; idx[4] = l; idx[5] = n; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = m; idx[4] = n; idx[5] = l; spin_batch.push_back( std::make_pair( idx, -val ) );

  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = m; idx[5] = n; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = n; idx[5] = m; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = n; idx[4] = m; idx[5] = l; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = n; idx[4] = l; idx[5] = m; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = m; idx[4] = l; idx[5] = n; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = m; idx[4] = n; idx[5] = l; spin_batch.push_back( std::make_pair( idx,  val ) );

  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = m; idx[5] = n; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = n; idx[5] = m; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = n; idx[4] = m; idx[5] = l; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = n; idx[4] = l; idx[5] = m; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = m; idx[4] = l; idx[5] = n; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = m; idx[4] = n; idx[5] = l; spin_batch.push_back( std::make_pair( idx, -val ) );

  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = m; idx[5] = n; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = n; idx[5] = m; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = n; idx[4] = m; idx[5] = l; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = n; idx[4] = l; idx[5] = m; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = m; idx[4] = l; idx[5] = n; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = m; idx[4] = n; idx[5] = l; spin_batch.push_back( std::make_pair( idx,  val ) );

  // Get transpose elements with same parity factors, hardcoded for speed
  if ( !skip_transpose ) {
    idx[0] = n; idx[1] = m; idx[2] = l; idx[3] = k; idx[4] = j; idx[5] = i; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx[0] = n; idx[1] = m; idx[2] = l; idx[3] = k; idx[4] = i; idx[5] = j; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = n; idx[1] = m; idx[2] = l; idx[3] = i; idx[4] = j; idx[5] = k; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = n; idx[1] = m; idx[2] = l; idx[3] = i; idx[4] = k; idx[5] = j; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx[0] = n; idx[1] = m; idx[2] = l; idx[3] = j; idx[4] = k; idx[5] = i; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = n; idx[1] = m; idx[2] = l; idx[3] = j; idx[4] = i; idx[5] = k; spin_batch.push_back( std::make_pair( idx,  val ) );
  
    idx[0] = n; idx[1] = l; idx[2] = m; idx[3] = k; idx[4] = j; idx[5] = i; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = n; idx[1] = l; idx[2] = m; idx[3] = k; idx[4] = i; idx[5] = j; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx[0] = n; idx[1] = l; idx[2] = m; idx[3] = i; idx[4] = j; idx[5] = k; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx[0] = n; idx[1] = l; idx[2] = m; idx[3] = i; idx[4] = k; idx[5] = j; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = n; idx[1] = l; idx[2] = m; idx[3] = j; idx[4] = k; idx[5] = i; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx[0] = n; idx[1] = l; idx[2] = m; idx[3] = j; idx[4] = i; idx[5] = k; spin_batch.push_back( std::make_pair( idx, -val ) );
  
    idx[0] = m; idx[1] = n; idx[2] = l; idx[3] = k; idx[4] = j; idx[5] = i; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = m; idx[1] = n; idx[2] = l; idx[3] = k; idx[4] = i; idx[5] = j; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx[0] = m; idx[1] = n; idx[2] = l; idx[3] = i; idx[4] = j; idx[5] = k; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx[0] = m; idx[1] = n; idx[2] = l; idx[3] = i; idx[4] = k; idx[5] = j; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = m; idx[1] = n; idx[2] = l; idx[3] = j; idx[4] = k; idx[5] = i; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx[0] = m; idx[1] = n; idx[2] = l; idx[3] = j; idx[4] = i; idx[5] = k; spin_batch.push_back( std::make_pair( idx, -val ) );
  
    idx[0] = m; idx[1] = l; idx[2] = n; idx[3] = k; idx[4] = j; idx[5] = i; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx[0] = m; idx[1] = l; idx[2] = n; idx[3] = k; idx[4] = i; idx[5] = j; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = m; idx[1] = l; idx[2] = n; idx[3] = i; idx[4] = j; idx[5] = k; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = m; idx[1] = l; idx[2] = n; idx[3] = i; idx[4] = k; idx[5] = j; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx[0] = m; idx[1] = l; idx[2] = n; idx[3] = j; idx[4] = k; idx[5] = i; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = m; idx[1] = l; idx[2] = n; idx[3] = j; idx[4] = i; idx[5] = k; spin_batch.push_back( std::make_pair( idx,  val ) );
  
    idx[0] = l; idx[1] = m; idx[2] = n; idx[3] = k; idx[4] = j; idx[5] = i; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = l; idx[1] = m; idx[2] = n; idx[3] = k; idx[4] = i; idx[5] = j; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx[0] = l; idx[1] = m; idx[2] = n; idx[3] = i; idx[4] = j; idx[5] = k; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx[0] = l; idx[1] = m; idx[2] = n; idx[3] = i; idx[4] = k; idx[5] = j; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = l; idx[1] = m; idx[2] = n; idx[3] = j; idx[4] = k; idx[5] = i; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx[0] = l; idx[1] = m; idx[2] = n; idx[3] = j; idx[4] = i; idx[5] = k; spin_batch.push_back( std::make_pair( idx, -val ) );
  
    idx[0] = l; idx[1] = n; idx[2] = m; idx[3] = k; idx[4] = j; idx[5] = i; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx[0] = l; idx[1] = n; idx[2] = m; idx[3] = k; idx[4] = i; idx[5] = j; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = l; idx[1] = n; idx[2] = m; idx[3] = i; idx[4] = j; idx[5] = k; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = l; idx[1] = n; idx[2] = m; idx[3] = i; idx[4] = k; idx[5] = j; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx[0] = l; idx[1] = n; idx[2] = m; idx[3] = j; idx[4] = k; idx[5] = i; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx[0] = l; idx[1] = n; idx[2] = m; idx[3] = j; idx[4] = i; idx[5] = k; spin_batch.push_back( std::make_pair( idx,  val ) );
  }

}

//===========================================================================================================================================================
// This is a general routine for arbitrary order permutations.
/*
void Fourpdm_permutations::get_even_and_odd_perms( const std::vector<int> mnpq, 
                                                   std::vector< std::vector<int> > & even_perms, 
                                                   std::vector< std::vector<int> > & odd_perms )
{
  // Get all even and odd mnpq permutations
  bool even = false;

  // Must sort them to get all possible permutations
  std::vector<int> foo = mnpq;
  std::sort( foo.begin(), foo.end() );

  // Get first set
  std::vector< std::vector<int> > perms1;
  do { 
    perms1.push_back( foo ); 
    if (foo == mnpq) even = true; 
  } while ( next_even_permutation(foo.begin(), foo.end()) );

  // Re-sort and swap LAST TWO elements to ensure we get all the remaining permutations
  std::sort( foo.begin(), foo.end() );
  assert( foo.size() == 4 );
  std::swap( foo[2], foo[3] );

  // Get second set
  std::vector< std::vector<int> > perms2;
  do { 
    perms2.push_back( foo ); 
  } while ( next_even_permutation(foo.begin(), foo.end()) );

  // Assign as even or odd permutations
  even_perms = perms1;
  odd_perms = perms2;
  if (!even) std::swap( even_perms, odd_perms );

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// This is a general routine for arbitrary order permutations.

void Fourpdm_permutations::get_spin_permutations_general( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
                                                          const std::vector<int>& indices, const double& val )
{
  assert( indices.size() == 8 );
  std::vector<int> idx;
  int i = indices[0];
  int j = indices[1];
  int k = indices[2];
  int l = indices[3];
  int m = indices[4];
  int n = indices[5];
  int p = indices[6];
  int q = indices[7];

  // If indices are not all unique, then all elements should be zero (and next_even_permutation fails)
  std::vector<int> v = {i,j,k,l};
  std::sort( v.begin(), v.end() );
  if ( (v[0]==v[1]) || (v[1]==v[2]) || (v[2]==v[3]) ) return;
  std::vector<int> w = {m,n,p,q};
  std::sort( w.begin(), w.end() );
  if ( (w[0]==w[1]) || (w[1]==w[2]) || (w[2]==w[3]) ) return;
  bool skip_transpose = ( v == w );

  // The number of possible combinations is (4!)**2 
  //------------------------------------------------

  // Get all even and odd permutations
  const std::vector<int> ijkl = {i,j,k,l};
  std::vector< std::vector<int> > ijkl_even; idx[] = ijkl_odd;
  get_even_and_odd_perms( ijkl; idx[] = ijkl_even, ijkl_odd );
  assert ( ijkl_even.size() + ijkl_odd.size() == 24 );

  const std::vector<int> mnpq = {m,n,p,q};
  std::vector< std::vector<int> > mnpq_even; idx[] = mnpq_odd;
  get_even_and_odd_perms( mnpq; idx[] = mnpq_even, mnpq_odd );
  assert ( mnpq_even.size() + mnpq_odd.size() == 24 );

  // Even-Even terms
  for ( auto u = ijkl_even.begin(); u != ijkl_even.end(); ++u ) {
    for ( auto v = mnpq_even.begin(); v != mnpq_even.end(); ++v ) {
      idx = { (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] };
      spin_batch.push_back( std::make_pair( idx, val ) );
      if ( !skip_transpose ) {
        // Include transpose
        std::reverse( idx.begin(); idx[] = idx.end() );
        spin_batch.push_back( std::make_pair( idx, val ) );
      }
    }
  }
  // Even-Odd terms
  for ( auto u = ijkl_even.begin(); u != ijkl_even.end(); ++u ) {
    for ( auto v = mnpq_odd.begin(); v != mnpq_odd.end(); ++v ) {
      idx = { (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] };
      spin_batch.push_back( std::make_pair( idx, -val ) );
      if ( !skip_transpose ) {
        // Include transpose
        std::reverse( idx.begin(); idx[] = idx.end() );
        spin_batch.push_back( std::make_pair( idx, -val ) );
      }
    }
  }
  // Odd-Even terms
  for ( auto u = ijkl_odd.begin(); u != ijkl_odd.end(); ++u ) {
    for ( auto v = mnpq_even.begin(); v != mnpq_even.end(); ++v ) {
      idx = { (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] }; 
      spin_batch.push_back( std::make_pair( idx, -val ) );
      if ( !skip_transpose ) {
        // Include transpose
        std::reverse( idx.begin(); idx[] = idx.end() );
        spin_batch.push_back( std::make_pair( idx, -val ) );
      }
    }
  }
  // Odd-Odd terms
  for ( auto u = ijkl_odd.begin(); u != ijkl_odd.end(); ++u ) {
    for ( auto v = mnpq_odd.begin(); v != mnpq_odd.end(); ++v ) {
      idx = { (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] };
      spin_batch.push_back( std::make_pair( idx, val ) );
      if ( !skip_transpose ) {
        // Include transpose
        std::reverse( idx.begin(); idx[] = idx.end() );
        spin_batch.push_back( std::make_pair( idx, val ) );
      }
    }
  }

}
*/

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//FIXME can speed up by replacing push_backs with [] ? Or .reserve() ?
// Note we can use the arbitary order routines above instead, but probably slower.

void Fourpdm_permutations::get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
                                                  const std::vector<int>& indices, const double& val )
{
  assert( indices.size() == 8 );
  std::vector<int> idx(8);
  int i = indices[0];
  int j = indices[1];
  int k = indices[2];
  int l = indices[3];
  int m = indices[4];
  int n = indices[5];
  int p = indices[6];
  int q = indices[7];

  // If indices are not all unique, then all elements should be zero (and next_even_permutation fails)
  std::vector<int> v(4);
  v[0] = i;
  v[1] = j;
  v[2] = k;
  v[3] = l;
  std::sort( v.begin(), v.end() );
  if ( (v[0]==v[1]) || (v[1]==v[2]) || (v[2]==v[3]) ) return;
  std::vector<int> w(4);
  w[0] = m;
  w[1] = n;
  w[2] = p;
  w[3] = q;
  std::sort( w.begin(), w.end() );
  if ( (w[0]==w[1]) || (w[1]==w[2]) || (w[2]==w[3]) ) return;
  bool skip_transpose = ( v == w );

  // The number of possible combinations is (4!)**2 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = l; idx[3] = j; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = j; idx[3] = k; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = l; idx[3] = k; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = i; idx[3] = l; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = k; idx[3] = i; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = j; idx[3] = l; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = l; idx[3] = i; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = i; idx[3] = j; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = k; idx[3] = j; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = i; idx[3] = k; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = j; idx[3] = i; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) );
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) );
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) );
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = m; idx[5] = n; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = m; idx[5] = p; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = m; idx[5] = q; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = n; idx[5] = m; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = n; idx[5] = p; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = n; idx[5] = q; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = p; idx[5] = m; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = p; idx[5] = n; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = p; idx[5] = q; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = q; idx[5] = m; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = q; idx[5] = n; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = q; idx[5] = p; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = j; idx[2] = l; idx[3] = k; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = k; idx[2] = j; idx[3] = l; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = i; idx[1] = l; idx[2] = k; idx[3] = j; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = i; idx[2] = k; idx[3] = l; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = k; idx[2] = l; idx[3] = i; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = j; idx[1] = l; idx[2] = i; idx[3] = k; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = i; idx[2] = l; idx[3] = j; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = j; idx[2] = i; idx[3] = l; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = k; idx[1] = l; idx[2] = j; idx[3] = i; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = i; idx[2] = j; idx[3] = k; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = j; idx[2] = k; idx[3] = i; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = m; idx[5] = n; idx[6] = q; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = m; idx[5] = p; idx[6] = n; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = m; idx[5] = q; idx[6] = p; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = n; idx[5] = m; idx[6] = p; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = n; idx[5] = p; idx[6] = q; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = n; idx[5] = q; idx[6] = m; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = p; idx[5] = m; idx[6] = q; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = p; idx[5] = n; idx[6] = m; idx[7] = q; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = p; idx[5] = q; idx[6] = n; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = q; idx[5] = m; idx[6] = n; idx[7] = p; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = q; idx[5] = n; idx[6] = p; idx[7] = m; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx[0] = l; idx[1] = k; idx[2] = i; idx[3] = j; idx[4] = q; idx[5] = p; idx[6] = m; idx[7] = n; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  // Get transpose elements with same parity factors, hardcoded for speed
  if ( !skip_transpose ) {
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = l; idx[5] = k; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = j; idx[5] = l; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = k; idx[5] = j; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = k; idx[5] = l; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = l; idx[5] = i; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = i; idx[5] = k; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = l; idx[5] = j; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = i; idx[5] = l; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = j; idx[5] = i; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = j; idx[5] = k; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = k; idx[5] = i; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = i; idx[5] = j; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = n; idx[3] = m; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = p; idx[3] = m; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = q; idx[3] = m; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = m; idx[3] = n; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = p; idx[3] = n; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = q; idx[3] = n; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = m; idx[3] = p; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = n; idx[3] = p; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = q; idx[3] = p; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = m; idx[3] = q; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = n; idx[3] = q; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = p; idx[3] = q; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = k; idx[5] = l; idx[6] = j; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = l; idx[5] = j; idx[6] = k; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = j; idx[5] = k; idx[6] = l; idx[7] = i; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = l; idx[5] = k; idx[6] = i; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = i; idx[5] = l; idx[6] = k; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = k; idx[5] = i; idx[6] = l; idx[7] = j; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = j; idx[5] = l; idx[6] = i; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = l; idx[5] = i; idx[6] = j; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = i; idx[5] = j; idx[6] = l; idx[7] = k; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = k; idx[5] = j; idx[6] = i; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = i; idx[5] = k; idx[6] = j; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = q; idx[2] = n; idx[3] = m; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = n; idx[2] = p; idx[3] = m; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = p; idx[2] = q; idx[3] = m; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = p; idx[2] = m; idx[3] = n; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = q; idx[2] = p; idx[3] = n; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = m; idx[2] = q; idx[3] = n; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = q; idx[2] = m; idx[3] = p; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = q; idx[1] = m; idx[2] = n; idx[3] = p; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = n; idx[2] = q; idx[3] = p; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = p; idx[1] = n; idx[2] = m; idx[3] = q; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = m; idx[1] = p; idx[2] = n; idx[3] = q; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx[0] = n; idx[1] = m; idx[2] = p; idx[3] = q; idx[4] = j; idx[5] = i; idx[6] = k; idx[7] = l; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  }
}

//===========================================================================================================================================================
//
//void Fourpdm_permutations::get_spatial_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
//                                                  const std::vector<int>& indices, const double& val )
//{
//  assert( indices.size() == 8 );
//  std::vector<int> idx;
//  int i = indices[0];
//  int j = indices[1];
//  int k = indices[2];
//  int l = indices[3];
//  int m = indices[4];
//  int n = indices[5];
//  int p = indices[6];
//  int q = indices[7];
//
//  // If indices are not all unique, then all elements should be zero (and next_even_permutation fails)
//  std::vector<int> v = {i,j,k,l};
//  std::sort( v.begin(), v.end() );
//  if ( (v[0]==v[1]) || (v[1]==v[2]) || (v[2]==v[3]) ) return;
//  std::vector<int> w = {m,n,p,q};
//  std::sort( w.begin(), w.end() );
//  if ( (w[0]==w[1]) || (w[1]==w[2]) || (w[2]==w[3]) ) return;
//  bool skip_transpose = ( v == w );
//
//  // The number of possible combinations is (4!)**2 
//  idx = { i, j, k, l, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
//
//}
//
//===========================================================================================================================================================

void Pairpdm_permutations::get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
                                                 const std::vector<int>& indices, const double& val )
{
  assert( indices.size() == 2 );
  std::vector<int> idx(2);
  int i = indices[0];
  int j = indices[1];

  idx[0] = i; idx[1] = j;
  spin_batch.push_back( std::make_pair( idx, val ) );
  // Transpose is same 
  if ( i != j ) {
    idx[0] = j; idx[1] = i;
    spin_batch.push_back( std::make_pair( idx, -val ) );
    // <m|a_ia_j|n>=-<m|a_ja_i|n>
  }

}
//===========================================================================================================================================================

}
