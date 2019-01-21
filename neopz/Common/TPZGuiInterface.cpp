/**
 * @file
 * @brief Contains the TPZGuiInterface methods.
 */
//$Id: TPZGuiInterface.cpp,v 1.3 2010-07-23 04:54:25 phil Exp $

#include "TPZGuiInterface.h"
#include "pzerror.h"
#include "TPZStream.h"
#include "Hash/TPZHash.h"

TPZGuiInterface::TPZGuiInterface(){
	this->fCanceled = false;
	//this->fMessage = "";
	this->fProgressBarPos = 0;
	this->fProgressBarMaxPos = 0;
	this->fProgressBarMinPos = 0;
}

TPZGuiInterface::~TPZGuiInterface(){
	//nothing to be done
}

int TPZGuiInterface::ClassId() const {
    return Hash("TPZGuiInterface");
}

void TPZGuiInterface::Read(TPZStream& buf, void* context) {
    buf.Read(fCanceled);
    buf.Read(&fMessage);
    buf.Read(&fProgressBarPos);
    buf.Read(&fProgressBarMaxPos);
    buf.Read(&fProgressBarMinPos);
}

void TPZGuiInterface::Write(TPZStream& buf, int withclassid) const {
    buf.Write(fCanceled);
    buf.Write(&fMessage);
    buf.Write(&fProgressBarPos);
    buf.Write(&fProgressBarMaxPos);
    buf.Write(&fProgressBarMinPos);
}

void TPZGuiInterface::UpdateCaption(){
    if(fMessage.length())
    {
        std::cout << fMessage.c_str() << "\n";
    }
    std::cout << "Progress bar = " << fProgressBarPos << "/" << fProgressBarMaxPos
	<< "\n";
}

/** Message attribute for UpdateGUI method */
std::string TPZGuiInterface::Message(){
    return fMessage;
}


/** Change the message of the object */
void TPZGuiInterface::SetMessage(const std::string &message)
{
    fMessage = message;
}


void TPZGuiInterface::Start(){
	std::cout << "Starting execution\n";
}

void TPZGuiInterface::End(){
	std::cout << "Execution finished\n";
}

void TPZGuiInterface::ShowErrorMessage(std::string message){
	PZError << message.c_str() << "\n";
//	DebugStop();
}


void TPZGuiInterface::SetKilled(){
	this->fCanceled = true;
}


bool TPZGuiInterface::AmIKilled(){
	return this->fCanceled;
}



