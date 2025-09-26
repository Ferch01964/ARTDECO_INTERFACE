def RaiseError(Log_file, Message):
    log_out = open(Log_file, "w")
    mesg = "\n\n###########\n##\n## CONFIGURATION ERROR\n##    |\n##    |_ "
    lines = Message.split("\n")
    for line in lines:
        mesg += line
        mesg += "\n##                      "
    mesg += "\n##\n###########\n"
    log_out.write(mesg)
    log_out.close()
    raise ValueError(mesg)
        
        
        
