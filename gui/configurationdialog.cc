#include "configurationdialog.h"
#include "cmbmainwindow.h"
#include "cmbmainwindow.h"

#include <QPixmap>
#include <QMessageBox>

#include <iostream>

#define PrinterProfileRole (Qt::UserRole+1)

static const char* const printerPixmatData[] = {
"64 64 551 2",
".c c None",
".# c None",
"Qt c None",
"fM c #000010",
"f1 c #000018",
"fN c #000020",
"dq c #000029",
"dy c #000031",
"fO c #000039",
"fT c #000041",
"gv c #00006a",
"gx c #000073",
"f9 c #00007b",
"aB c #0024ff",
"aj c #002cff",
"aC c #0044ff",
"fP c #080839",
"ai c #0838c5",
"aU c #0838ff",
"aV c #0861e6",
"f0 c #101041",
"fH c #10184a",
"aA c #1038ff",
"#4 c #103cff",
"aT c #1040c5",
"az c #1044c5",
"b. c #1048cd",
"#5 c #1079ff",
"aD c #1091ff",
"gj c #181441",
"dP c #181808",
"gD c #181c52",
"f8 c #18205a",
"ah c #1848ff",
"#P c #184cff",
"#Q c #185dd5",
"b# c #1895ff",
"fB c #20245a",
"#O c #2050cd",
"#3 c #2050ff",
"ay c #2055ff",
"bt c #205dff",
"ak c #208dff",
"dX c #292418",
"cu c #292c62",
"#2 c #2959d5",
"ag c #2959ff",
"af c #295dcd",
"bs c #295dd5",
"a9 c #295dff",
"ax c #2961d5",
"#1 c #297dff",
"bu c #2981ee",
"e. c #313029",
"ae c #3161ff",
"aS c #3165d5",
"bH c #316dff",
"aw c #3185ff",
"#R c #3195ff",
"em c #393831",
"aQ c #396dd5",
"aR c #396dff",
"br c #3975de",
"bI c #39b6ff",
"el c #413c39",
"eF c #41406a",
"es c #41446a",
"d0 c #414473",
"do c #41447b",
"fm c #414483",
"dS c #4144bd",
"gp c #41487b",
"ge c #41488b",
"dd c #4148bd",
"gE c #4148c5",
"dk c #414c7b",
"b0 c #414c83",
"aY c #414c8b",
"dt c #414cc5",
"d7 c #414ccd",
"e6 c #41508b",
"bF c #4175ff",
"a6 c #4179ff",
"bq c #417dff",
"bV c #4181ff",
"a7 c #4199ff",
"eA c #4a4441",
"c1 c #4a4c7b",
"de c #4a507b",
"cT c #4a5083",
"by c #4a508b",
"cQ c #4a50c5",
"eT c #4a50cd",
"fd c #4a50d5",
"bc c #4a558b",
"#T c #4a5594",
"aE c #4a5994",
"fL c #4a599c",
"bp c #4a81e6",
"bU c #4a85ff",
"b9 c #4a91ff",
"bW c #4aa1f6",
"eL c #524c4a",
"cS c #52558b",
"cR c #525983",
"cH c #52598b",
"gh c #525994",
"#x c #52599c",
"gc c #5259d5",
"#7 c #525d94",
"fc c #525d9c",
"eN c #525dac",
"gy c #525dde",
"gz c #526194",
"#G c #52619c",
"gm c #5261a4",
"bo c #5289ff",
"bn c #528de6",
"a4 c #528dff",
"cn c #5291ee",
"b7 c #5291ff",
"av c #5295ff",
"#6 c #52baff",
"eX c #5a5552",
"eW c #5a5952",
"ct c #5a5dd5",
"du c #5a618b",
"ca c #5a6194",
"am c #5a61de",
"gJ c #5a6594",
"eQ c #5a659c",
"#m c #5a65a4",
"gn c #5a65de",
"#y c #5a65e6",
"al c #5a699c",
"go c #5a69a4",
"#t c #5a6dbd",
"eq c #5a71c5",
"fk c #5a71ff",
".6 c #5a75cd",
"bm c #5a95ff",
"bT c #5a9dee",
"aP c #5a9dff",
"co c #5aa1ff",
"#0 c #5aaaff",
"aW c #5adaff",
"dw c #5aff08",
"e9 c #62615a",
"dQ c #62656a",
"dM c #626994",
"gu c #62699c",
"gr c #6269a4",
"#U c #6269e6",
"eH c #626da4",
"#z c #626dac",
"#k c #6271ac",
"fb c #6271b4",
"eE c #6271f6",
"eG c #6275b4",
".E c #6275cd",
"e1 c #6279c5",
".f c #6279cd",
"ec c #627dcd",
"eb c #627dd5",
"cm c #629dff",
"bB c #62a1ff",
"cl c #62a5f6",
"bS c #62a5ff",
"ad c #62b2ff",
"dv c #62b66a",
"#M c #62baff",
"fg c #6a696a",
"ff c #6a6d6a",
"d8 c #6a6d9c",
"gF c #6a71a4",
"ev c #6a71ac",
"g. c #6a71ee",
"e4 c #6a75b4",
"ee c #6a75bd",
"#w c #6a75ee",
"eO c #6a75ff",
"an c #6a79b4",
"dZ c #6a79c5",
"aG c #6a79ff",
"dA c #6a7dbd",
"cd c #6a7dc5",
"fZ c #6a7dcd",
".S c #6a7dd5",
".2 c #6a7dff",
"fY c #6a81cd",
".h c #6a81d5",
".0 c #6a81ff",
"fl c #6a85d5",
"bR c #6aaaf6",
"cC c #6aaaff",
"cN c #6aaef6",
"cE c #6ab2ff",
"#Z c #6abeff",
"ba c #6ac2ff",
"bQ c #6aceff",
"#N c #6ae6ff",
"dW c #732820",
"fq c #737173",
"ey c #7375a4",
"eU c #737dac",
"ei c #737db4",
"#V c #737dbd",
"#o c #737dff",
"#q c #7381c5",
"#I c #7381ff",
"cx c #7385c5",
".a c #7385ff",
"fF c #7389d5",
".N c #7389de",
"fS c #7389ff",
"#C c #739dff",
"b5 c #73b2ff",
"bl c #73b6f6",
"cj c #73b6ff",
"b4 c #73baff",
"aO c #73beff",
"ab c #73c2ff",
"ac c #73e6ff",
"gH c #7b3cb4",
"gi c #7b40b4",
"fG c #7b40bd",
"fx c #7b797b",
"fC c #7b7d7b",
"dO c #7b7d8b",
"cs c #7b81b4",
"d3 c #7b81ff",
"eB c #7b85a4",
"e# c #7b85ac",
"fe c #7b85b4",
"d4 c #7b85bd",
"#g c #7b85cd",
"fu c #7b85ff",
"cJ c #7b89c5",
"b2 c #7b89cd",
"a8 c #7b89ff",
".l c #7b8dd5",
"bG c #7b8dff",
".Y c #7b91d5",
"fA c #7b91de",
".e c #7b91ff",
"cA c #7bbaff",
"cy c #7bc2ff",
"#X c #7bc6ff",
"#D c #7bcaff",
"#B c #7bceff",
"ci c #7bdeff",
"bv c #7beaff",
"ga c #837d73",
"f5 c #83817b",
"fw c #838183",
"dV c #83858b",
"eY c #8389a4",
"#i c #838dd5",
"#s c #838dff",
".t c #8391cd",
".A c #8391ff",
"eD c #8395d5",
"bM c #8395de",
"d1 c #8395e6",
"eo c #8395ff",
"bP c #83c6ff",
"as c #83caff",
"c7 c #83ceff",
"db c #83deee",
"cM c #83eaff",
"a# c #83eeff",
"et c #8b4cbd",
"cI c #8b4cc5",
"dF c #8b50c5",
"fW c #8b898b",
"gb c #8b8d9c",
"dT c #8b8dcd",
"eS c #8b95bd",
"fa c #8b95d5",
"eu c #8b95ff",
"gd c #8b99d5",
"dG c #8b99e6",
"gC c #8b99ff",
"c4 c #8b9de6",
".d c #8b9dee",
"a5 c #8b9dff",
"bE c #8ba1ff",
"b8 c #8ba5ff",
"bJ c #8bceff",
"aq c #8bd2ff",
"dp c #9475c5",
"c3 c #9475cd",
"bL c #9479cd",
"#F c #947dde",
"fh c #9495a4",
"bK c #9495ff",
"ex c #9499bd",
"bN c #9499cd",
"dI c #9499d5",
"gB c #94a1de",
".n c #94a1ff",
"cb c #94a5f6",
"cv c #94aaff",
"aL c #94d6ff",
"eK c #9c7594",
"gK c #9c7dcd",
"c2 c #9c7dd5",
"dl c #9c81d5",
"gt c #9c81de",
"aF c #9c85de",
"#8 c #9c85e6",
"ek c #9c959c",
"dE c #9c99b4",
"fQ c #9c9d9c",
"d6 c #9c9dbd",
"dL c #9ca1bd",
"d. c #9ca1c5",
"cP c #9ca1cd",
".H c #9ca5d5",
"dB c #9ca5ff",
".r c #9caaee",
"fE c #9caaff",
"bD c #9cb2ff",
"c# c #9cd2ff",
"a2 c #9cdaff",
"dD c #a485de",
"bd c #a485e6",
"gf c #a489e6",
"#l c #a489ee",
"dY c #a48dd5",
"#H c #a48dee",
"#n c #a48df6",
"fs c #a491ff",
"#a c #a499ff",
"#v c #a4aacd",
"dr c #a4aaee",
"gG c #a4aaff",
"be c #a4aecd",
"#e c #a4b2e6",
"bC c #a4beff",
"cD c #a4c2ff",
"c. c #a4d2ff",
"c8 c #a4daff",
"a0 c #a4deff",
"e8 c #ac85ac",
"dU c #ac8de6",
"f2 c #ac95f6",
".m c #ac99ff",
".R c #ac9dff",
"eJ c #acaaac",
"fy c #acaeb4",
"gA c #acb2d5",
"c5 c #acb2ee",
"dh c #acb2ff",
"gw c #acb6d5",
"fi c #acb6de",
"f# c #acb6ff",
".j c #acbeff",
"b6 c #acceff",
"cF c #acd2ff",
"bi c #ace2ff",
"dm c #acea94",
"fV c #b475b4",
"#A c #b479ff",
".g c #b481ff",
"e2 c #b485ff",
"#f c #b489ff",
"f4 c #b48dac",
"ej c #b495e6",
"f3 c #b499f6",
"b3 c #b49dff",
".F c #b4a1ff",
".p c #b4a5ff",
"#E c #b4b6ff",
"gq c #b4bad5",
".O c #b4baff",
"gI c #b4bee6",
".y c #b4beff",
"c6 c #b4c6ff",
"da c #b4ceee",
"dj c #b4ceff",
"ck c #b4d2ff",
"cB c #b4d6ff",
"a. c #b4daff",
"bg c #b4deff",
"#Y c #b4e2ff",
"#L c #b4e6ff",
"fU c #bd85f6",
"#h c #bd89ff",
"cc c #bd8dff",
"ea c #bd91ff",
"fp c #bd99bd",
"eI c #bd9dee",
"dg c #bda1ff",
"cG c #bda5f6",
"#p c #bda5ff",
"f7 c #bdaaff",
"f6 c #bdaeff",
"ep c #bdb2ff",
"g# c #bdbab4",
"fJ c #bdbabd",
"cK c #bdbaff",
"#9 c #bdbede",
"dH c #bdbef6",
"cW c #bdbeff",
"gk c #bdc2d5",
".4 c #bdc2ff",
"fz c #bdc6d5",
"ce c #bdc6e6",
"dR c #bdc6f6",
"ao c #bdc6ff",
"cL c #bdceff",
"a3 c #bddaff",
"#K c #bddeff",
"cZ c #bde2ff",
"aa c #bde6ff",
"au c #bdeeff",
".s c #c591ff",
"e7 c #c5a5f6",
"eM c #c5aae6",
"en c #c5aaee",
"eh c #c5aaff",
"cw c #c5aeff",
".b c #c5b6ff",
"di c #c5c2ff",
"#b c #c5c6d5",
"fj c #c5c6e6",
"ds c #c5c6ff",
".8 c #c5cad5",
".T c #c5cade",
"ed c #c5caf6",
"#W c #c5caff",
"aH c #c5cee6",
"fR c #c5ceff",
"cz c #c5e6ff",
"aX c #c5eaff",
"at c #c5eeff",
"cg c #c5f2ff",
"fn c #cdaeff",
"e3 c #cdb6ff",
"b1 c #cdbaff",
"dz c #cdbeff",
"ft c #cdc2ff",
"fo c #cdcacd",
".Z c #cdcad5",
"dJ c #cdcaff",
"## c #cdceff",
"fX c #cdd2de",
".w c #cdd2ee",
".G c #cdd2f6",
"d2 c #cdd2ff",
".P c #cdd6ff",
"aZ c #cddaff",
"cO c #cddeff",
"#u c #cde2f6",
"cp c #cde2ff",
"bY c #cde6ff",
"bw c #cdeaff",
"bO c #cdeeff",
"ar c #cdf2ff",
"aN c #cdf6ff",
"d9 c #d58dde",
"fI c #d599ff",
"eC c #d5a1ff",
".i c #d5a5ff",
"gg c #d5b2ee",
"f. c #d5b6e6",
"e5 c #d5b6ff",
"bb c #d5baff",
"df c #d5beff",
"cU c #d5c2ff",
".k c #d5c6ff",
".3 c #d5d2d5",
"d5 c #d5d2ff",
"#. c #d5d6e6",
"#r c #d5d6ee",
"dc c #d5d6ff",
"c9 c #d5deff",
"aM c #d5f2ff",
"cf c #d5faff",
"gs c #debeff",
".z c #dec2ff",
".7 c #dec6ff",
"#d c #dedaee",
"dK c #dedaff",
"fK c #dedede",
".X c #dedeee",
"eR c #dedeff",
"#J c #dee2ee",
".J c #dee2f6",
"bf c #dee2ff",
"bX c #def2ff",
"bk c #def6ff",
"aK c #defaff",
"ez c #e6a1e6",
"fr c #e6c2ee",
"c0 c #e6c6ff",
"fv c #e6caff",
"eZ c #e6ceff",
"er c #e6d2ff",
"#c c #e6deee",
".1 c #e6e2ee",
"eP c #e6e2ff",
".M c #e6e6ee",
"aJ c #e6e6f6",
"cr c #e6e6ff",
"dN c #e6eaff",
"ch c #e6faff",
"a1 c #e6ffff",
"#S c #eeceff",
".I c #eed2ff",
"gl c #eed6ff",
".q c #eedaff",
"e0 c #eee6f6",
"dx c #eee6ff",
".L c #eeeaf6",
"cq c #eeeaff",
".K c #eeeef6",
"bZ c #eeeeff",
".9 c #eef2ff",
"bj c #eeffff",
"fD c #f6d2f6",
"eV c #f6d2ff",
"bz c #f6d6ff",
".U c #f6daff",
"aI c #f6eef6",
"ew c #f6eeff",
".W c #f6f2f6",
"bx c #f6f2ff",
".5 c #f6f6f6",
".V c #f6f6ff",
"ap c #f6faff",
"bh c #f6ffff",
"cX c #ffbaff",
"cV c #ffc2ff",
"dC c #ffc6ff",
"cY c #ffcaff",
"ef c #ffceff",
"bA c #ffd2ff",
"eg c #ffd6ff",
"d# c #ffdaff",
".B c #ffdeff",
".D c #ffe2ff",
".C c #ffe6ff",
"#j c #ffeaff",
".Q c #ffeeff",
".u c #fff2ff",
".v c #fff6ff",
".x c #fffaff",
"dn c #ffff6a",
".o c #ffffff",
"Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.a.b.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.d.e.f.g.h.i.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.j.k.l.m.n.o.o.p.aQt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.c.q.r.s.t.u.o.o.v.o.w.g.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#.x.y.z.A.B.v.o.o.o.C.o.D.o.E.F.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.cQt.cQt.cQt.G.C.H.I.J.o.o.o.K.o.L.o.M.o.o.I.NQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.#Qt.D.v.O.I.P.o.o.o.v.o.u.o.Q.o.C.o.D.o.o.R.SQt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.c.o.T.U.T.o.o.o.o.o.V.o.W.o.K.o.L.o.M.o.X.o.Y.a.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#.o.Q.v.Z.o.o.o.o.o.x.o.x.o.v.o.u.o.Q.o.C.o.C.o.u.o.0Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.x.o.x.o.o.o.o.o.x.o.x.o.v.o.W.o.K.o.L.o.M.o.1.o.o.2.SQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#.o.3.o.o.o.o.o.o.o.o.o.x.o.v.o.u.o.Q.o.C.o.C.o.B.o.4.p.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.c.x.5.o.o.o.o.o.o.o.x.o.V.o.W.o.K.o.L.o.M.o.1.o.1.o.6.0.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.D.v.o.o.o.o.o.o.o.o.x.o.v.o.u.o.Q.o.C.o.C.o.B.o.o.7.aQt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.c.o.8.o.o.o.o.o.o.o.x.o.v.o.W.o.K.o.L.o.M.o.1.o.X.o.9.2.EQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#.u.D.o.o.o.o.o.o.o.x.o.v.o.u.o.Q.o.C.o.C.o.B.o.B.o.a.p.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt#..C.o.o.o.o.o.o.x.o.V.o.W.o.K.o.L.o.M.o.1.o.X.o.v.v.hQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt##.o.o.o.o.o.o.o.x.o.v.o.u.o.Q.o.C.o.C.o.B.o.B.o.o#a.0Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.c.v#b.o.o.o.o.o.x.o.v.o.W.o.K.o.L.o.M.o.1.o#c.o#d.o#e#f.cQt#g#h#iQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.u.Q.x.o.o.o.o.o.x.o.v.o.u.o.Q.o.C.o.C.o.B.o.D.o#j.o.f.F#k#l#m#n#o#p#qQt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.cQt#r.v.o.o.x.o.x.o.V.o.W.o.K.o.L.o.M.o.1.o.X.o.X.o.o#s#t.B#u.o#v#w#x#y#z#A.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.#.x.O.o.o.o.o.o.x.o.v.o.u.o.Q.o.C.o.C.o.B.o.B.o.B.o#j.o#B#C#D.o.o.o#E#F#G#H#IQt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.cQt.c.C#J.o.x.o.x.o.v.o.W.o.K.o.L.o.M.o.1.o#c.o.1.o#K#L#M#N#O#P#Q#R#L.o.o#S#T#U#VQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.#Qt#W.B.o.o.o.o.x.o.v.o.u.o.Q.o.C.o.C.o.B.o.C.o#X#Y#Z#L#0#1#2#3#4#5#6.o.o.o#7#8#oQt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.cQt.cQt#9.o.o.o.x.o.V.o.W.o.K.o.L.o.M.o.1.oa.a##Xaaabacadaeafag#Oahaiajak.o.o.oalamanQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.#Qt.#.Qao.o.x.o.x.o.v.o.u.o.Q.o#j.o.Dapaqarasat#Dauavawax.2ay#3azaAaBaCaD.o.o.oaEaFaGQt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.cQt.cQtaH.I.o.o.x.o.v.o.W.oaI.oaJ.oa.aKaLaMaqaNaOaPaQaRaSae#2#3aTaUaVaWaX.o.o.o.o.oaYam.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.#Qt.#Qt.O.v.o.o.x.o.v.o.v.oaZ.oa0a1a2aKaLa3a4a5a6a7aRa8a9#3b.b#ba.o.x.o.o.o.v.o.obbbcbd.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.cQt.cQt.c.Dbe.o.o.o.x.obf.vbgbhbibja0bkblbmbnbobpbqbraRbsbtbubvbw.o.o.o.o.o.V.obx.obxby#7Qt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#bz#j.o.o.obA.oau.o#Lbha2a3bBbCbmbDbobEbFbGbHbIbJ.o.x.o.o.o.o.o.x.o.v.o.Q.obKbLbcQt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.cQt.cQtbM#IbNbzbf.ubO.oat.o#LbhbPbQbRbSbTbmbnbUbrbVbWbXbY.o.o.o.o.o.o.o.x.o.V.obZ.obZ.ob0am.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.#b1b2b3#I.u.P.o.x.oat.o#ja1b4a3b5b6bSbCb7b8b9c.c#.o.o.o.o.o.o.o.o.o.x.o.v.o.u.o#j.o.xbLcaQt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQtcbcccd#oce#ja1.ocf.ocg.oaXchabcicjckclcmcnco#Z.9cp.o.o.o.o.o.o.o.o.o.x.o.V.o.9.ocq.ocr.ocsctcuQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qtcvcwcx.B.o.o.o.o.D.oaN.o.QbjcyczcAcBcCcDcEaXcF.o.o.o.o.o.o.o.o.o.o.o.o.o.x.o.u.o.Q.o#j.o.vcGcHcI.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.ccccJ.DcKbz#E.o.o.ocL.oat.obPcMcycBcNcibJbxcO.o.o.o.o.o.o.o.o.o.o.o.o.o.x.obx.obZ.o.9.ocPcQcRcScTQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#cUcxbzcV.DcW.DcX.C.o.ocY.oaucZ#XcZasaKck.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.x.o.v.o.u.o.vc0c1c2cK.ob0c3.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQtc4.ac5.CcW.D.4.DcW.C.o.oc6.oc7bkc8.Qc9.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.x.o.V.o.x.od.b0cHd#dadbdcdddeQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#df#Idgdh#jcV.Cdi.CcV.u.o.ocVbZdj.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.x.o.x.o.xc0dkdlcX.o#jdmdn.ododpdqQt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQtc4.gc5aGdr.Qds.Cds#j.4.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.od.dtdu.D.V.ocr.odvdwdxdoc1dy.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#dzdA.u#E.mdB.ucY#j#W#jdC.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.oc0b0dDds.o.u.o.u.odE.ocO.QdodFdqQt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQtdG.2#W.QdHaGdI.udJ.Q###jdK.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.odLb0dM.Q.o.obZ.odNdOdPdQdN.odRdSc1dy.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#dzaG.QdJ.ucYb3dT.xbA.u###j.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.o.oc0b0dUbA.o.x.o.x.odVdWdXdYdZ.o.D#jd0dFdqQt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQtd1aG#W.QdJ.ud2d3d4.vd5.u##.o.o.o.o.o.o.o.o.o.o.o.o.o.o.od6d7d8.x.o.o.V.o.9d9e.e.e#eaebaGec.uedd0c1dy.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#b1ee.Qef.Q##.vegehei.xeg.v##.o.o.o.o.o.o.o.o.o.o.o.o.zaYejd#.o.o.o.o.oekelemeneoep.P.o.ob1eqeresetdqQt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQtbMaG##.u##.ud5.vdKeuev.udK.vew.o.o.o.o.o.o.o.o.oexaYey.o.o.o.o.oapezeAeAeBeCeD.x.o.oaJ.o.o.U.6eEeFQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#b1eG.ud5.ubA.vd5.x.D.zeH.Q.B.o.o.o.o.o.o.o.xbbaYeI.D.o.o.o.o.oeJeKeLeMdBcU.B.o.u.o.C.o.B.o.o.oeq.peNQt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQtbMeOd2.vd5.vdc.vdc.xeP.IeQ#jeR.o.o.o.o.oeSeTeU.o.o.o.o.o.oeVeWeXeYeZ.H.o.o.o.K.oe0.o.M.o.X.o.o.oe1e2bEQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#e3e4.D.B.xdc.xeg.xd#.o.C.D#G.B.o.o.ve5e6e7.C.o.o.o.o.ocXe8e9f.f#.I.C.o.v.o.Q.o#j.o.C.o.D.o.B.o.o.ofa.R.eQt.#Qt.#Qt.#Qt.#",
".cQt.cQt.d#Ifb.BeP.xdK.xdK.oeR.ocr.Qfcd#eSfdfe.o.o.o.o.o.o.Dfffgfh.Dfi.o.o.o.V.o.W.o.K.o.L.oaJ.o.X.o.X.o.K.ofjfkflQt.cQt.cQt.cQt",
"Qt.#Qt.#Qtfm#pfbbz.C.od#.o.B.o.B.o#j.xeQfn#j.o.o.o.o.ofofpfqfrds.B#j.o.o.o.v.o.u.o.Q.o#j.o.C.o.D.o.B.o.B.o.D.o.ufs.aft.#Qt.#Qt.#",
".cQt.cQt.cQtcufu#zfvcq.oeR.obf.oeP.ocr.o.o.o.o.o.o.Qfwfxfy.ufz.o.o.o.x.o.x.o.V.o.W.o.K.o.L.o.M.o#c.o.X.o.X.o.V.o.o.2fAQt.cQt.cQt",
"Qt.#Qt.#Qt.#QtfBdg#z.7.Q.o.D.o.D.o.C.o.Q.o.o.o.o.BfCfDd#.Q.u.o.o.o.x.o.x.o.v.o.u.o.Q.o#j.o.C.o.D.o.B.o.B.o.o.ofE.RfFfG.#Qt.#Qt.#",
".cQt.cQt.cQt.cQtfH#I#zfIbZ.ocr.ocr.ocr.o.o.o.o.ofJ.ofK.o.o.o.o.o.o.o.x.o.x.o.V.o.W.o.K.o.L.oaJ.o.X.o.L.o.o.U.feafLfMfNfO.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.#QtfP.F#zfn.u.o.C.o#j.o.u.o.o.o.o.BfQ.od#.o.o.o.o.o.x.o.x.o.v.o.u.o.Q.o#j.o.C.o.B.o.o.ofR.RfS.mfNfNdyfT.#Qt.#Qt.#",
".cQt.cQt.cQt.cQt.cQtfOaG#zfU.V.ocq.ocq.o.o.o.o.ofJfVfW.ofX.o.o.o.o.o.x.o.x.o.V.o.W.o.K.oe0.oaJ.o.o.ofYfSfZf0f1dydyQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.#Qt.#Qtdyf2#zf3.v.o.Q.o.v.o.o.o.od#fqf4f5.o##.o.o.o.x.o.x.o.v.o.u.o.Q.o#j.o.o.o.uf6.hf7f8f1dqf9.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.cQt.cQt.cQt.cg.eHg..V.obZ.o.o.o.o.og#gagbgc#7.o#b.o.o.o.x.o.x.o.V.o.K.o.L.o.o.ogd.2fFgefMdqdyfT.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.##HeHgf.u.o.x.o.o.o.oggghaFeQgigj.ogk.o.o.o.x.o.v.o.u.o.x.o.ogl.0f6gmf1dqfOfOQt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.cQt.cQt.cQt.cQt.cgngoambZ.o.o.oeSaY#7#UgpfNfNdqdy.ogq.o.o.o.v.o.W.o.o.ocea8bMeOfMfNdqfO.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#dlgr#F#jgsaYgtgucIf1gvdyfO.#Qt.#.ogw.o.o.o.x.o.o.o.tdfa8dyf1gxfOQt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.ccSgrgygz#UdkfNfNdydyQt.cQt.cQt.c.ogA.x.o.o.LeCgBgCgDf1dqfO.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#gEgFdFf1dqdyf9.#Qt.#Qt.#Qt.#Qt.#.odh.ugGeZfEgHf1dqdyQt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cdqfNdydyQt.cQt.cQt.cQt.cQt.cQt.c.ugI.qgJfMfNfOfOQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt",
"Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#fO.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#gKfMgvdyfT.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#Qt.#",
".cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cdydyQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt.cQt"};

ConfigurationDialog::ConfigurationDialog(vector<LowLevelPlot>& profiles, QWidget* parent , const char* name , bool modal)  : QDialog(parent), Ui::DesignConfigurationDialog() , PrtProfiles(profiles) , Point(2.835) {

  setupUi(this);

//X   listBox->setIconSize(QSize(listBox->width()-12, listBox->width()-12));
//X   QListWidgetItem * item = new QListWidgetItem(listBox);
//X   QPixmap printerPixmap(printerPixmatData);
//X   QIcon printerIcon(printerPixmap.scaledToWidth(listBox->width()));
//X   item->setIcon(printerIcon);

  CurrentProfile = -1;

  XAxisMarks->addItem("decimal");
  YAxisMarks->addItem("decimal");
  XAxisMarks->addItem("exponent");
  YAxisMarks->addItem("exponent");

  connect( saveButton, SIGNAL( clicked() ), this, SIGNAL( saveProfile() ) );
  connect( listBox, SIGNAL( itemActivated( QListWidgetItem* ) ), this, SLOT( activateProfile( QListWidgetItem* ) ) );

  listBox->setContextMenuPolicy( Qt::ActionsContextMenu );
  QAction * deleteProfile = new QAction( "Delete Profile", listBox );
  listBox->addAction( deleteProfile );
  connect( deleteProfile, SIGNAL( triggered() ), this, SLOT( deleteProfile() ) );
}

void ConfigurationDialog::deleteProfile()
{
  if ( listBox->currentItem()->data( PrinterProfileRole ).toInt() < 7 ) {
    QMessageBox::information( this, "Information", "You cannot delete the standard printing profiles from this list.",
		   QMessageBox::Ok );
    return;
  }
  PrtProfiles.erase( PrtProfiles.begin()+listBox->currentItem()->data( PrinterProfileRole ).toInt() );
  updateProfiles();
}

void ConfigurationDialog::updateProfiles() {
  listBox->clear();
  cout << "PrtProfiles to update: " << PrtProfiles.size() << endl;
  for ( unsigned int i = 0; i < PrtProfiles.size(); ++i ) {
    cout << "adding profile: " << i << endl;
    QListWidgetItem * item = new QListWidgetItem(listBox);
    item->setData( PrinterProfileRole, i );
    item->setText( QString::fromStdString( PrtProfiles[ i ].Name ) );
  }
}

void ConfigurationDialog::activateProfile( QListWidgetItem* item ) {
  if ( !item ) return;
  setProfile( item->data(PrinterProfileRole).toInt() );
}

void ConfigurationDialog::setProfile(int idx) {
  if (idx > (int)PrtProfiles.size() -1 ) throw Bad_Error("configurationdialog:: no such plot"); 
  cout << "CONFIGURATIONDIALOG:: Setting profile to: " << idx << endl;
  syncProfile();
  CurrentProfile  = idx;
  LowLevelPlot &p = PrtProfiles[idx];
  safety = p; // make a copy
  cout << "Now synchronizing"<<endl;
  syncDialog(p);
  if ( listBox->count() == 0 ) return;
  for ( int i=0; i <= listBox->count(); ++i ) {
    if ( listBox->item( i ) && ( listBox->item( i )->data( PrinterProfileRole ).toInt() == idx) )
      listBox->setItemSelected( listBox->item( i ), true );
  }
}

void ConfigurationDialog::syncDialog(LowLevelPlot& p) {
  XPaper->setText(CmbMainWindow::toStr(p.Width,1));
  YPaper->setText(CmbMainWindow::toStr(p.Height,1));
  AxisLabelSize->setValue(p.AxisLabelSize);
  TickLabelSize->setValue(p.TickLabelSize);
  XOffset->setText(CmbMainWindow::toStr(p.XLabelOffset,4));
  YOffset->setText(CmbMainWindow::toStr(p.YLabelOffset,4));
  AutomaticTick->setChecked(p.AutomaticTick);
  XAxisMarks->setCurrentIndex(labelStyle2int(p.XTickLabelStyle));
  YAxisMarks->setCurrentIndex(labelStyle2int(p.XTickLabelStyle));

  LeftMargin->setValue((int)rint(p.LeftMargin*100));
  RightMargin->setValue((int)rint(p.RightMargin*100));
  TopMargin->setValue((int)rint(p.TopMargin*100));
  BottomMargin->setValue((int)rint(p.BottomMargin*100));

  FrameLineWidth->setValue(p.FrameLineWidth);
  CurveLineWidth->setValue(p.CurveLineWidth);
  TicksLineWidth->setValue(p.TicksLineWidth);

  XLabelText->setText(p.XLabelText.c_str());
  YLabelText->setText(p.YLabelText.c_str());

  XPrecision->setValue(p.Significant_x);
  YPrecision->setValue(p.Significant_y);

  QString aString;
  StartX->setText(aString.setNum(p.StartTick_x));
  StartY->setText(aString.setNum(p.StartTick_y));
  XSpacing->setText(aString.setNum(p.StepTick_x));

  cout <<"config compare: ";
  cout.precision(16);
  cout << p.StepTick_x << " " << p.float2string(p.StepTick_x) << endl;

  YSpacing->setText(aString.setNum(p.StepTick_y));

  setWindowTitle(p.Name.c_str());

}

void ConfigurationDialog::syncProfile() {
  if (CurrentProfile > (int)PrtProfiles.size() -1 || CurrentProfile < 0) return ; 
  cout << "First, saving old profile: " << CurrentProfile << endl;

  LowLevelPlot &p = PrtProfiles[CurrentProfile];

  p.Width = XPaper->text().toFloat();
  p.Height = YPaper->text().toFloat();
  p.AxisLabelSize = AxisLabelSize->value();
  p.TickLabelSize = TickLabelSize->value();
  p.XLabelOffset = XOffset->text().toFloat();
  p.YLabelOffset = YOffset->text().toFloat();
  p.AutomaticTick = AutomaticTick->isChecked();
  p.StartTick_x = StartX->text().toFloat();
  p.StartTick_y = StartY->text().toFloat();
  p.StepTick_x = XSpacing->text().toFloat();
  p.StepTick_y = YSpacing->text().toFloat();

  p.XTickLabelStyle = int2labelStyle(XAxisMarks->currentIndex());
  p.YTickLabelStyle = int2labelStyle(YAxisMarks->currentIndex());

  p.LeftMargin=LeftMargin->text().toFloat()*0.01;
  p.RightMargin=RightMargin->text().toFloat()*0.01;
  p.TopMargin=TopMargin->text().toFloat()*0.01; 
  p.BottomMargin=BottomMargin->text().toFloat()*0.01;

  p.FrameLineWidth = FrameLineWidth->value();
  p.CurveLineWidth = CurveLineWidth->value();
  p.TicksLineWidth = TicksLineWidth->value();

  p.XLabelText = XLabelText->text().toStdString();
  p.YLabelText = YLabelText->text().toStdString();

  p.Significant_x = XPrecision->value();
  p.Significant_y = YPrecision->value();
}

void ConfigurationDialog::reject() {
  cout << "I REJECT" << endl;
  PrtProfiles[CurrentProfile] = safety;
  syncDialog(safety);
  QDialog::reject();
}

void ConfigurationDialog::accept() {
  syncProfile();
  QDialog::accept();
}

LowLevelPlot::LabelStyle ConfigurationDialog::int2labelStyle(int i) {
  switch (i) {
  case 0:
    return LowLevelPlot::decimal;
  default:
    return LowLevelPlot::exponent;
  }
}

int ConfigurationDialog::labelStyle2int(LowLevelPlot::LabelStyle s) {
  switch (s) {
  case LowLevelPlot::decimal:
    return 0;
  default:
    return 1;
  }
}

