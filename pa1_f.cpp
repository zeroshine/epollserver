#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <stack>
#include <sys/time.h>
#define Minsup 300
using namespace std;
class node{
    public:
	map<string,node*> children;
	node* parent;
	node* mapnode;
	string item;
	int freq;
	int index;
	int headindex;
	vector<node*> link;
	vector<string> lis;


};
class Comp{
    public:
	bool operator ()(const node* ln,const node* rn)	{
	    return ln->freq>rn->freq;
	}

};
class Comps{
    public:
	bool operator ()(const node* ln,const node* rn)	{
	    return ln->item<rn->item;
	}

};
void buildtree(vector<vector<node*> >table,vector<vector<node*> > &lisvec,vector<string> vn){
    map<string,node*> first_map;
    node* null=new node;
    for(vector<vector<node*> >::iterator vvit=table.begin();vvit!=table.end();++vvit){
	node* lastnode = null;
	for(vector<node*>::iterator vit=(*vvit).begin();vit!=(*vvit).end();++vit){
	    node* n = *vit;
	    map<string,node*>::iterator mfit=lastnode->children.find(n->item);
	    map<string,node*>::iterator firstit=first_map.find(n->item);
	    if(mfit==lastnode->children.end()){
		node* treenode=new node;
		treenode->item=n->item;
		treenode->parent=lastnode;
		treenode->index=n->index;
		treenode->headindex=n->index;
		lastnode->children.insert(pair<string,node*>(n->item,treenode));
		if(firstit==first_map.end()){
		    first_map.insert(pair<string,node*>(n->item,treenode));
		    treenode->link.push_back(treenode);
		}else{
		    firstit->second->headindex+=n->index;
		    firstit->second->link.push_back(treenode);
		}
		lastnode=treenode;
	    }else{
		mfit->second->index+=n->index;
		firstit->second->headindex+=n->index;
		lastnode=mfit->second;
	    }
	}
    }
    for(map<string,node*>::iterator mit=first_map.begin();mit!=first_map.end();++mit){
	vector<node*> linkvec=mit->second->link;
	vector<vector<node*> > newtable;
	for(vector<node*>::iterator vit=linkvec.begin();vit!=linkvec.end();++vit){
	    if((*vit)->parent!=null){
		node* tempnode=(*vit)->parent;
		vector<node*> v;
		int tempindex=(*vit)->index;
		stack<node*> st;
		while(tempnode!=null){
		    node* storenode=new node;
		    storenode->item=tempnode->item;
		    storenode->index=tempindex;
		    st.push(storenode);
		    tempnode=tempnode->parent;
		}

		while(!st.empty()){
		    v.push_back(st.top());
		    st.pop();
		}
		newtable.push_back(v);
	    }
	}
	node* lisn=new node;
	for(vector<string>::iterator vit=vn.begin();vit!=vn.end();++vit){
	    lisn->lis.push_back(*vit);
	}
	if(mit->second->headindex>=Minsup){
	    lisn->lis.push_back(mit->first);
	    lisvec[lisn->lis.size()-1].push_back(lisn);
	}
	if(newtable.size()!=0){
	    buildtree(newtable,lisvec,lisn->lis);
	}
    }
}
int main(){
    timeval timev,pretimev;
    gettimeofday(&pretimev,NULL);
    ifstream infile("hw1data/1000.txt");
    string line;
    map<string ,node*> first_map;
    vector< vector<node*> > table;
    vector< vector<node* > >lisvec;
    while(getline(infile,line)){
	string word;
	istringstream ss(line);
	vector<node*> v;
	while(ss>>word&&word!="!EOD"){
	    node* tnode =new node;
	    v.push_back(tnode);
	    tnode->index=0;
	    tnode->freq=1;
	    tnode->item=word;
	    map<string,node*>::iterator fmit=first_map.find(word);
	    if(fmit==first_map.end()){
		node* newnode=new node;
		newnode->item=word;
		newnode->freq=1;
		newnode->index=0;
		tnode->mapnode=newnode;
		first_map.insert(pair<string,node*>(word,newnode));
	    }else{
		tnode->mapnode=fmit->second;
		++(fmit->second->freq);
	    }
	}
	sort(v.begin(),v.end(),Comps());
	table.push_back(v);
    }
    for(vector<vector<node*> >::iterator vvit=table.begin();vvit!=table.end();++vvit){
	for(vector<node*>::iterator vit=(*vvit).begin();vit!=(*vvit).end();){
	    (*vit)->freq=(*vit)->mapnode->freq;
	    if((*vit)->freq<Minsup){
		(*vvit).erase(vit);
	    }else{
		++vit;
	    }
	}
	stable_sort((*vvit).begin(),(*vvit).end(),Comp());
    }
    for(map<string,node*>::iterator mit=first_map.begin();mit!=first_map.end();){
	node* n = mit->second;
	if(n->freq<Minsup){
	    delete n;
	    first_map.erase(mit++);
	}else{
	    ++mit;
	}
    }
    lisvec.reserve(first_map.size());
    for(int i=0;i<first_map.size();++i){
	vector<node*> v;
	lisvec.push_back(v);
    }
    node* null=new node;
    null->item="null";
    for(vector<vector<node*> >::iterator vvit=table.begin();vvit!=table.end();++vvit){
	node* lastnode = null;
	for(vector<node*>::iterator vit=(*vvit).begin();vit!=(*vvit).end();++vit){
	    node* n = *vit;
	    map<string,node*>::iterator mfit=lastnode->children.find(n->item);
	    if(mfit==lastnode->children.end()){
		node* treenode=new node;
		treenode->item=n->item;
		treenode->parent=lastnode;
		treenode->index=1;
		lastnode->children.insert(pair<string,node*>(n->item,treenode));
		(first_map.find(n->item))->second->link.push_back(treenode);
		lastnode=treenode;
	    }else{
		++(mfit->second->index);
		lastnode=mfit->second;
	    }
	}
    }
    for(map<string,node*>::iterator mit=first_map.begin();mit!=first_map.end();++mit){
	vector<node*> linkvec=mit->second->link;
	vector<vector<node*> > newtable;
	newtable.reserve(linkvec.size());
	for(vector<node*>::iterator vit=linkvec.begin();vit!=linkvec.end();++vit){
	    if((*vit)->parent!=null){
		vector<node*> v;
		node* tempnode=(*vit)->parent;
		int tempindex=(*vit)->index;
		stack<node*> st;
		while(tempnode!=null){
		    node* storenode=new node;
		    storenode->item=tempnode->item;
		    storenode->index=tempindex;
		    st.push(storenode);
		    tempnode=tempnode->parent;
		}
		while(!st.empty()){
		    v.push_back(st.top());
		    st.pop();
		}
		newtable.push_back(v);
	    }
	}
	vector<string> vn;
	vn.push_back(mit->first);
	buildtree(newtable,lisvec,vn);
    }
    gettimeofday(&timev,NULL);
    float diff=(timev.tv_sec*1000000+timev.tv_usec-pretimev.tv_sec*1000000-pretimev.tv_usec)/1000000.000000;
    cout<<"total time: "<<diff<<endl;
    for(int i=0;i<lisvec.size();++i){
	cout<<"large item set "<<i+1<<" : "<<lisvec[i].size()<<endl;
    }
}
