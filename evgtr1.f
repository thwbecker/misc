c
c
c     subroutine to determine tectonic regions based on jordan's (1981) GTR1 regionalization
c
c
c     obtained on 08/15/06 from Bogdan Kustowski
c
c     $Id: gtr1_code.f,v 1.1 2006/08/15 17:06:47 becker Exp becker $
c
c
c
      subroutine evgtr1(xlat,xlon,itype)
      character*72 gtr1(36)
      character*1 type
      data gtr1(1)  /'ccccccccccccccccccccccbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbcccccccc'/
      data gtr1(2)  /'cccccccccccccccqqqqqqqqqqqqqqqqqqqbbabbbbbbqqqqqqqqqqqbbbbbbbbbbcccccccc'/
      data gtr1(3)  /'cccccccccccqqqqqqqqqqqsssssssssqqqbbaaqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq'/
      data gtr1(4)  /'qqqqccccccqqpppppppssbbbqsssssqqbbbaaabqqqqqqqqqqqqqpppppppppqqqqqqqqqqq'/
      data gtr1(5)  /'qqqqqqqqqpppsssssssssssbbssssqbbaabbbqqqssssppqqqppppppppspppqqqqqqqqqqq'/
      data gtr1(6)  /'qqqqqqqqqqqppsssssppsssqbbssbbaabbbbqqsssssppppqpppppppppppppppqqqqqqqqq'/
      data gtr1(7)  /'cqqqqqbbqqqpppssspppsssqbbbbaabbbbqqqqsppppppppqppppppspppqqqqqqqqqqqccc'/
      data gtr1(8)  /'qqqbbbbbaaqqqpppssspsssssqbbbabbbqqqqqpppppppppqqqqppqqqqqqqqqqqqqqqbbqq'/
      data gtr1(9)  /'bbbbbbbbaaaqqqpppsssssqqqqqbbbabbbqqqqqqqpspppppqqqqqqqqqqqqppqqqbqccbbb'/
      data gtr1(10) /'bbbbbbbbbaaqqqqppppppqqccccbbaabbbqqqbqqqqbbqqppqqqqqqqqqqqqpqbbqqcccccb'/
      data gtr1(11) /'bbbbbbbbbbbqqqqpppqqqcccccbbbaabbcqqbbqcqqqqqqqqqqqqpqqqqpppqqbqqccccccc'/
      data gtr1(12) /'cccbbbbbbbbbqqqqqqqqqcccbbbaabbbccqqpppcpccppqqqqqpqqqqqqqqpqqqbqccccccc'/
      data gtr1(13) /'cccbbbbbbbbbaaqqqccqqccbbbaabbbccqppppppppqpppqqqqpppqqqqppqqqbbqccccccc'/
      data gtr1(14) /'ccbbbbbbbbbbaaaqqcqqqccbbbaabbbcqsppsspspppasppqbqppsqqqqpqqqbbbqccccccc'/
      data gtr1(15) /'ccccbbbbbbbbaaaqqqqaqccqbbaabbbcqqppppppsppsappbbbqsqcbqpqbbqbbbbqcccccc'/
      data gtr1(16) /'ccccbbbbbbbbaaaabbqqcccqbbaabbbcqqpsssppsssqaaaabbbscbbqqqbbqbbbqccccccc'/
      data gtr1(17) /'ccccbbbbbbbbaaaaabbqqqpsqqbbaabbbqsssssssspqppbbabbbcbbqqqbqqqbbbbcccccc'/
      data gtr1(18) /'ccccbbbbbbbaaaaaaaaaqpssssbbbaabbbbbbcsppssqpbbbaabbbbbqqqqqbqbbbbcccccc'/
      data gtr1(19) /'ccccbbbbbbbaaaaaaabbqppppppqbbbbaabbbcpppqsqcbbbbabbbbbbqqqqqqqqqbqccccc'/
      data gtr1(20) /'ccccbbbbbbbaaaaaaabbqpppssppscbbbaabbcqpsssscbbbbabbbbbbbqqqqqqqqqbqcccc'/
      data gtr1(21) /'ccccbbbbbbbaaaaaaabbqqpsssssccbbaabbbcqspsssccbbaabbbbbbbbbcqqqqqbbbbqaa'/
      data gtr1(22) /'aqccbbbbbbbaaaaaaabbbqqpppsscbbbaaabbcqppssqcqbbaabbbbbbbbcqpspsqqbbbqaa'/
      data gtr1(23) /'qccbbbbbbbaaaaaaaabbbbqpppsqcbbbaabbbbcppspccqbbaaabbbbbbcqsspssqqqbbbbb'/
      data gtr1(24) /'qcccbbbbbbaaaaaaaabbbbqpppqcbbbbaabbbbcppsqccbbbaaaabbbbbcqsspppqqqbbbbb'/
      data gtr1(25) /'qccccbbbbbaaaaaaaaabbqqppqcbbbbbaabbbbccqqccbbbabbaaabbbbccssqqsqqqbbbbb'/
      data gtr1(26) /'cccccbbbbbbaaaaaaaaabqppqccbbbbbaabbbbbcccbbbaabbbaaaaabbbbbbbbbqqbbbbbq'/
      data gtr1(27) /'ccccbbbbbbbaaaaabaaaaqqqqcbbbbbbaabbbbbbbbbaabbbbbbaaaaaaabaaabbbqbbbqqb'/
      data gtr1(28) /'cbbbbbbbbbbaaaaabbaaaqqqcccbbbbbaaabbbbbbbabbbbbbbbbaaaaaaaaaaaabbbbbqqq'/
      data gtr1(29) /'bbbbbbbabbaaaabbbbbbaqqqqqqbbbbbbbaaaaaaaabbbbbbbbbbbbbaaaaaaaaaabbbbqqb'/
      data gtr1(30) /'bbbbaaaaaabbbbbbbbbbbaaaaaabbaqbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbaaabbbbb'/
      data gtr1(31) /'baaaaabbbbbbbbbbbbbbbbaqqqccccbbbbbbbbbbbbbccccccccccccbbbbbbbbbbbaaaaab'/
      data gtr1(32) /'aabbbbbbbbbbbbbbbbbbbqqqqcccccccccccccccccccqsssqqqqsssppppppppppqqqbbaa'/
      data gtr1(33) /'bbbbbbbbbbbbqqqqqqqqqqqqqqqcccccqqssssssssssssssssssppspppppppssssqqqqqq'/
      data gtr1(34) /'qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqsspppppsspppppppppppssssppppppppppqqqqqq'/
      data gtr1(35) /'qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqssssssssssssssssssssspppppppppppqqqqqqqq'/
      data gtr1(36) /'qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqppppppppppppppppppppppppppppppppppqqqqqqqq'/


      ilat=1+int((90.-xlat)/5.)
      ilon=1+int((xlon+360.)/5.)
      ilon=ilon-36
      if(ilat.eq.37) ilat=36
      if(ilon.gt.72) ilon=ilon-72
      if(ilat.lt.1.or.ilat.gt.36) then
        write(6,"('latitude out of range',f8.2,i4)") xlat,ilat
        itype = -2
        return
      endif
      if(ilon.lt.1.or.ilon.gt.72) then
        write(6,"('longitude out of range',f8.2,i4)") xlon,ilon
        itype = -3
        return
      endif
c
      type=gtr1(ilat)(ilon:ilon)
      if(type.eq.'a') then
        itype=1
      else if(type.eq.'b') then
        itype=2
      else if(type.eq.'c') then
        itype=3
      else if(type.eq.'p') then
        itype=4
      else if(type.eq.'q') then
        itype=5
      else if(type.eq.'s') then
        itype=6
      else 
         print *, 'error in evgtr1'
         itype = -1
      endif
      return
      end
