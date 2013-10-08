#include <iostream>
#include "add.h"
#include <sys/types.h>          /* See NOTES */
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/epoll.h>
#include <errno.h>
#include <map>
#include <vector>
#define MAX_EVENTS 1000
using namespace std;
void setnonblocking(int sock){
    int opts;
    opts=fcntl(sock,F_GETFL);
    if(opts<0){
	cout<<"fcntl F_GETFL"<<endl;
    }
    opts = opts|O_NONBLOCK;
    if(fcntl(sock,F_SETFL,opts)<0){
	cout<<"fcntl F_SETFL"<<endl;
    } 
}
int main()
{
    ////////////////////client connect to Server/////////////////////////////////////////////////////
    int sockfd_tr;
    struct sockaddr_in dest_tr;
    char buffer_tr[1024];
    //char resp[10]="clientack";

    /* create socket */
    sockfd_tr = socket(PF_INET, SOCK_STREAM, 0);

    /* initialize value in dest */
    bzero(&dest_tr, sizeof(dest_tr));
    dest_tr.sin_family = PF_INET;
    dest_tr.sin_port = htons(8889);
    dest_tr.sin_addr.s_addr = inet_addr("140.112.94.58");
    //inet_aton("127.0.0.1", &dest.sin_addr);

    /* Connecting to server */
    int con=connect(sockfd_tr, (struct sockaddr*)&dest_tr, sizeof(dest_tr));

    /* Receive message from the server and print to screen */
    cout<<"connect to server "<<con<<endl;
    bzero(buffer_tr, 1024);
    //cout<<recv(sockfd_tr,buffer_tr,122,0)<<buffer_tr<<endl;
    ////////////////////Server act as PDC///////////////////////////////////
    struct sockaddr_in dest;
    char recieve[275];
    int clientfd,insockfd;
    struct sockaddr_in client_addr;
    socklen_t addrlen = sizeof(client_addr);
    int epfd,nfds;
    struct epoll_event ev,events[MAX_EVENTS];
    epfd = epoll_create(MAX_EVENTS);

    /* create sockett */
    int sockfd = socket(PF_INET, SOCK_STREAM, 0);
    setnonblocking(sockfd);
    setnonblocking(sockfd_tr);
    /* initialize structure dest */
    bzero(&dest, sizeof(dest));
    dest.sin_family = AF_INET;
    dest.sin_port = htons(9000);
    /* this line is different from client */
    dest.sin_addr.s_addr = INADDR_ANY;

    /* Assign a port number to socket */
    bind(sockfd, (struct sockaddr*)&dest, sizeof(dest));

    /* make it listen to socket with max 20 connections */
    listen(sockfd, MAX_EVENTS);

    ev.data.fd=sockfd;
    ev.events=EPOLLIN|EPOLLET;
    int s=epoll_ctl(epfd,EPOLL_CTL_ADD,sockfd,&ev);
    if(s==-1){
	cout<<"err sockfd"<<endl;
    }
    //int j =1;
    vector<map <long long int,Alg*>* > v;
    v.reserve(100);
    for(int i=0;i<100;i++){
	map<long long int,Alg*>* m=new map<long long int,Alg*>;
	v.push_back(m);
    }
    while(1){
	nfds=epoll_wait(epfd,events,MAX_EVENTS,-1);
	cout<<"Wait success "<<nfds<<endl;
	for(int i=0;i<nfds;i++){
	    cout<<"nfds : "<<nfds<<" i "<<i<<endl;
	    cout<<"now fd is"<<events[i].data.fd<<endl;
	    if ((events[i].events & EPOLLERR) ||(events[i].events & EPOLLHUP) ||
		    (!(events[i].events & EPOLLIN)))
	    {
		cout<<"epoll error"<<endl;
		close (events[i].data.fd);
		continue;
	    }

	    if(events[i].data.fd==sockfd){
		cout<<"catch some connect fd"<<events[i].data.fd<<endl;
		while(1){
		    cout<<"Wait accept"<<endl;
		    clientfd = accept(sockfd, (struct sockaddr*)&client_addr, &addrlen);
		    if(clientfd==-1){
			if ((errno == EAGAIN)||(errno == EWOULDBLOCK))
			{
			    break;
			}else{
			    break;
			}
		    }
		    setnonblocking(clientfd);
		    cout<<"connect from "<<inet_ntoa(client_addr.sin_addr)<<" fd "<<clientfd<<endl;
		    ev.data.fd=clientfd;
		    ev.events=EPOLLIN|EPOLLET;
		    if(epoll_ctl(epfd,EPOLL_CTL_ADD,clientfd,&ev)==-1){
			cout<<"epoll add error"<<endl;
		    }
		}
		cout<<"accept new down"<<endl;
		continue;	
	    }else if(events[i].events & EPOLLIN) {
		cout<<events[i].data.fd<<"fd can be read "<<endl;
		int read_byte;
		bool isR=false;
		getpeername(events[i].data.fd,(struct sockaddr*)&client_addr,&addrlen);
		string ip_from =inet_ntoa(client_addr.sin_addr);
		int port_from =ntohs(client_addr.sin_port);
		cout<<"socket from "<<ip_from<<" port: "<<port_from<<endl;
		if(ip_from=="140.112.145.218"){
		    cout<<"=====================L A B 1======================"<<endl;
		}else if(ip_from=="140.112.145.219"){
		    cout<<"=====================L A B 2======================"<<endl;
		}else if(ip_from=="140.112.145.220"){
		    cout<<"=====================L A B 3======================"<<endl;
		}else if(ip_from=="140.112.145.221"){
		    cout<<"=====================L A B 4======================"<<endl;
		}else{
		    cout<<"@@@@@@@@@@@@@@ N B @@@@@@@@@@@@@@"<<endl;
		    isR=true;
		}
		if((insockfd=events[i].data.fd)<0) continue;
		cout<<"insockfd > 0"<<endl;
		clientfd=events[i].data.fd;
		string temp;
		while(1){
		    bzero(recieve,275);
		    read_byte = recv(clientfd,recieve,1,0);
		    if(read_byte<0){
			if(errno!=EAGAIN){
			    cout<<"errno"<<endl;
			    cout<<"err:"<<strerror(errno)<<endl;
			    events[i].data.fd=-1;
			    break;
			}
			break;
			/*
			   close(insockfd);
			   events[i].data.fd=-1;
			   break;
			   */
		    }else if(read_byte==0){
			cout<<"close socket from "<<inet_ntoa(client_addr.sin_addr)<<endl;
			close(insockfd);
			events[i].data.fd=-1;
			break;
                    }else if(strcmp(recieve,"E")!=0){
			//cout<<recieve<<endl;
                        temp.append(recieve);
		    }else if(strcmp(recieve,"E")==0){
			string rcopy = temp;
			bzero(recieve,275);
			temp.clear();
			cout<<"read_byte: "<<rcopy.length()<<endl;
			//			cout<<"read time: "<<j<<endl<<endl;
			cout<<"recieve message: "<<rcopy<<endl<<endl;
			Alg* com=new Alg;
		       //	cout<<"claim object Alg"<<endl;
			com->parse(rcopy,isR);
			//cout<<"parse string and location_index is "<<com->getlocation_index()<<endl;
			map<long long int,Alg*>* datamap=v[com->getlocation_index()];
			//cout<<"claim map"<<endl;
			map<long long int,Alg*>::iterator it=datamap->find(com->gettimestamp());
			//cout<<"claim iterator"<<endl;
			if(it==datamap->end()){
			  //  cout<<"insert data"<<endl;
			    datamap->insert(pair<long long int,Alg*>(com->gettimestamp(),com));
			}else{
			    delete com;
			    it->second->parse(rcopy,isR);
			   // cout<<"parse down"<<endl;
			    it->second->compute();
			    cout<<"###################"<<it->second->getType()<<"########"<<endl;
			    //cout<<"compute transmission line parameters"<<endl<<endl;

			    it++;
			    /*
			    for(map<long long int,Alg*>::iterator its=datamap->begin();its!=it;++its){
				if(its!=datamap->end()){
				    delete its->second;
				    its->second=0;
				    datamap->erase(its);
				}else{
			//	    cout<<"delete end"<<endl;
				}
			    }
			    */
			}
			//string forsend = com.compute();
			//cout<<"compute transmission line parameters"<<endl<<endl;
			//char s[1024];
			//int send_len=strlen(forsend.c_str());
			//cout<<"data len"<<send_len<<endl;
			//sprintf(s,"%d %s",send_len,forsend.c_str());
			//cout<<"send byte: "<<send(sockfd_tr,forsend.c_str(),strlen(forsend.c_str()),0)<<endl;
			//cout<<"send message:"<<forsend<<endl<<endl;
			/*
			   char s[3];
			   sprintf(s,"%d",j);
			   cout<<"send ack "<<write(insockfd,s,6)<<endl;
			   */
			//			j++;
		    }
		}
	    }
	}
    }
}
