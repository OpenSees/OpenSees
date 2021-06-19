#/void\* OPS_.*()/s|.*|#ifdef OPS_API_COMMANDLINE\n&|;tx; b; :x; n; s|^\}|}\n#endif|; Tx
#/^void \*OPS_[A-z]*()/s|.*|#ifdef OPS_API_COMMANDLINE\n&|;tx; b; :x; n; s|^\}|}\n#endif|; Tx
#/^void\* OPS_[A-z]*(void)/s|.*|#ifdef OPS_API_COMMANDLINE\n&|;tx; b; :x; n; s|^\}|}\n#endif|; Tx
#/^void \*OPS_[A-z]*(void)/s|.*|#ifdef OPS_API_COMMANDLINE\n&|;tx; b; :x; n; s|^\}|}\n#endif|; Tx
#/^OPS_[A-z]*(void)/s|.*|#ifdef OPS_API_COMMANDLINE\n&|;tx; b; :x; n; s|^\}|}\n#endif|; Tx
/^OPS_[A-z]*()/s|.*|#ifdef OPS_API_COMMANDLINE\n&|; tx; b; :x; n; s|^\}|}\n#endif|; Tx
