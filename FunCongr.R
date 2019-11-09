library("ape")
library("tree")
library("gtools")
library("phangorn") #Packages needed for the analysis

isEmpty <- function(x) { #This function checks if a data frame is empty or not
  return(length(x)==0)
}

options(expressions=1e5) #This option is necessary to perform a great number of permutations

tree=read.tree("gtdb_r86.ssu.bacteria.fasttree_name_pruneado_bueno.tree") #The tree is loaded here
tab_info_pair=read.table("tab_info_pair.txt") #Within this table we have all the distances values (Jaccard and 16S) for each  posible pair.
tab_info_pair$comb<-paste(tab_info_pair$V1,tab_info_pair$V2,sep = "")

####With this function we are able to get all the information from the nodes and leaves, their ID and given name

node_labels_in_edge <- tree$node.label[tree$edge[,1]-Ntip(tree)]
tips_nodes <- tree$edge[,2]

select.tip.or.node <- function(element, tree) {
  ifelse(element < Ntip(tree)+1, tree$tip.label[element], tree$node.label[element-Ntip(tree)])
}

edge_table <- data.frame(
  "parent" = tree$edge[,1],
  "par.name" = sapply(tree$edge[,1], select.tip.or.node, tree = tree),
  "child" = tree$edge[,2],
  "chi.name" = sapply(tree$edge[,2], select.tip.or.node, tree = tree)
)

edge_table$parent=as.numeric(edge_table$parent)
edge_table$child=as.numeric(edge_table$child)
edge_table$par.name=as.character(edge_table$par.name)
edge_table$chi.name=as.character(edge_table$chi.name)

####

nodes_list=unique(edge_table$par.name) #We get all possible node names from the table

tab_random=data.frame(matrix(ncol=1001)) 

tab_random=rbind(tab_random,c(5,unlist(sapply(seq_len(1000),function(x) mean(sample(tab_info_pair$V4,5)))))) #We originate a initial vector with random values for 5 leaves. This vector will be used to assess the ecdf function

tab_random=tab_random[-1,] #We erase the blank row line

bad_leaves=NULL #Within this object we will save the tips from each node of the tree that we considere as bad nodes, the nodes that didn't pass the ecdf function check

means_node_tab=data.frame(matrix(ncol=3)) #In this table we will save the info of the node name, the mean of the node and the total number of pairwise combinations
p=1

for(i in 1:length(nodes_list)) #This loop is the core of the script. It will go over all nodes form last to first one, saving the info of the means and saving the leaves of the "bad nodes"
{
  leaves=Descendants(tree,edge_table[edge_table$par.name==nodes_list[length(nodes_list)-i+1],1])[[1]] #We save the information of all leaves of a certain node
  descendant_nodes=Descendants(tree,edge_table[edge_table$par.name==nodes_list[length(nodes_list)-i+1],1],type = "child")[[1]] #Knowing which are the children of each node we can get the info of means and number of pairwise for each descendant indepently
  descendant_nodes1=Descendants(tree,descendant_nodes[1])[[1]]
  descendant_nodes2=Descendants(tree,descendant_nodes[2])[[1]]
  
  if(length(leaves)>=2 & length(na.omit(match(leaves,bad_leaves)))==0) #Here we check if a certain node is multifurcated
  {
    print(paste(Sys.time(),nodes_list[length(nodes_list)-i+1]))
    
    if(isEmpty(tab_random[tab_random$X1==length(leaves),2])) #If a random distribution doesn't exist for a certain number of combinations, we create that distribution
    {
      tab_random=rbind(tab_random,c(length(leaves),unlist(sapply(seq_len(1000),function(x) mean(sample(tab_info_pair$V4,length(leaves)))))))
    }
    
    val_ecdf=ecdf(tab_random[tab_random$X1==length(leaves),2:1001]) #The ecdf funcion is created here
    
    vector1=NULL #Here we save the leaves names from each descendant
    vector2=NULL
    
    for(o in 1:length(descendant_nodes1))
    {
      vector1=c(vector1,edge_table[edge_table$child==descendant_nodes1[o],4])  
    }
    
    for(o in 1:length(descendant_nodes2))
    {
      vector2=c(vector2,edge_table[edge_table$child==descendant_nodes2[o],4])  
    }
    
    vector1=sort(vector1) #Leaves names are sorted here
    vector2=sort(vector2)
    
    leaves_combinations=data.frame(expand.grid(vector1,vector2)) #We perform the combinations and save them as characters
    leaves_combinations$Var1=as.character(leaves_combinations$Var1)
    leaves_combinations$Var2=as.character(leaves_combinations$Var2)
    
    values_j=NULL #Within this vector we will save all jaccard values from the pairwise of all leaves
  
    combinations_parse<-paste(leaves_combinations$Var1,leaves_combinations$Var2,sep = "") #Here we extract all the jaccard values from al the combinationes of leaves calculated
    values_j<-tab_info_pair[tab_info_pair$comb %in% combinations_parse,4]
    
    combinations_parse<-paste(leaves_combinations$Var2,leaves_combinations$Var1,sep = "")
    values_j<-c(values_j,tab_info_pair[tab_info_pair$comb %in% combinations_parse,4])

    if(length(leaves)>2 & (length(descendant_nodes1)>1 | length(descendant_nodes2)>1)) #Here we create a vector of values from the means and number of pairwise we got from all descendants and the new combinations created beetween each descendant
    {
      values_j=as.vector(c(values_j,rep(means_node_tab[means_node_tab$X1==descendant_nodes[1],2],means_node_tab[means_node_tab$X1==descendant_nodes[1],3]),rep(means_node_tab[means_node_tab$X1==descendant_nodes[2],2],means_node_tab[means_node_tab$X1==descendant_nodes[2],3])))
    }
    
    if(val_ecdf(mean(values_j))>0.05 & length(leaves)>=5) #Here we check if our node passes the threeshold of the ecdf function. In case the p-value obtained is greater of 0.05, the node is considerated as "bad" and its leaves are saved to avoid the evaluation of the subsequents nodes that posseses the sames leaves
    {
      print(paste(length(nodes_list)[1]-i+1,nodes_list[length(nodes_list)[1]-i+1],mean(values_j),sd(values_j))) #In the output, we will be able to identify a "bad node" as they will show the mean and sd of the computed jaccard
      bad_leaves=c(bad_leaves,leaves)
    }
    
    means_node_tab[p,1]=edge_table[edge_table$par.name==nodes_list[length(nodes_list)-i+1],1][1] #For each node we save their name
    means_node_tab[p,2]=mean(values_j) #And the mean of their jaccard values
    if(length(leaves)>2 & (length(descendant_nodes1)>1 | length(descendant_nodes2)>1)) #Here we save its number of pairwise, checking if the node has one or two descendants, adding the number of pairs of internal combinations to the number of pairs found within each child node
    { 
      if(!isEmpty(means_node_tab[means_node_tab$X1==descendant_nodes[1],3])&!isEmpty(means_node_tab[means_node_tab$X1==descendant_nodes[2],3]))
      {
        means_node_tab[p,3]=dim(leaves_combinations)[1]+means_node_tab[means_node_tab$X1==descendant_nodes[2],3]+means_node_tab[means_node_tab$X1==descendant_nodes[1],3]
      }else if(!isEmpty(means_node_tab[means_node_tab$X1==descendant_nodes[1],3]))
      {
        means_node_tab[p,3]=dim(leaves_combinations)[1]+means_node_tab[means_node_tab$X1==descendant_nodes[1],3]
      }else if(!isEmpty(means_node_tab[means_node_tab$X1==descendant_nodes[2],3]))
      {
        means_node_tab[p,3]=dim(leaves_combinations)[1]+means_node_tab[means_node_tab$X1==descendant_nodes[2],3]
      }
    }else
    {
      means_node_tab[p,3]=dim(leaves_combinations)[1] #If the node is multifurcated, we only keep the number of pairwise generated in that node
    }
    
    p=p+1
  }
}
