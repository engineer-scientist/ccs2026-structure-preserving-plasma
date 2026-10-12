mkdir -p ~/grokbot-sp/p3 && cd ~/grokbot-sp/p3 && tar xzf /tmp/gb_p3.tgz && rm /tmp/gb_p3.tgz
for f in ex2 ex2_326; do [ -x ../$f ] && ln -sf ../$f $f; done
ls -la; ldd ./ex2* 2>/dev/null | grep -i petsc
