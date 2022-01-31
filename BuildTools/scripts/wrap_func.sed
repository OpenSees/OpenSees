# Claudio Perez
# this sed script wraps a C function in an ifdef block
#           1
# 01234567890123456789
# /<LINE_PATTERN>/s|<FIND_PATTERN>|<REPLACE_PATTERN>\n&|
# for LINE in file:
#     if match(LINE, LINE_PATTERN);

#/void\* OPS_.*()/s|.*|#ifdef OPS_API_COMMANDLINE\n&|;tx; b; :x; n; s|^\}|}\n#endif|; Tx
#/^void \*OPS_[A-z]*()/s|.*|#ifdef OPS_API_COMMANDLINE\n&|;tx; b; :x; n; s|^\}|}\n#endif|; Tx
#/^void\* OPS_[A-z]*(void)/s|.*|#ifdef OPS_API_COMMANDLINE\n&|;tx; b; :x; n; s|^\}|}\n#endif|; Tx
#/^void \*OPS_[A-z]*(void)/s|.*|#ifdef OPS_API_COMMANDLINE\n&|;tx; b; :x; n; s|^\}|}\n#endif|; Tx
#/^OPS_[A-z]*(void)/s|.*|#ifdef OPS_API_COMMANDLINE\n&|;tx; b; :x; n; s|^\}|}\n#endif|; Tx
/^OPS_[A-z]*()/s|.*|#ifdef OPS_API_COMMANDLINE\n&|; tx; b; :x; n; s|^\}|}\n#endif|; Tx
