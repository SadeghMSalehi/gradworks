//
//  QGraphicsImageItem.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/7/13.
//
//

#include "QGraphicsImageItem.h"


unsigned int __grays[256] = {
    4278190080u, 4278255873u, 4278321666u, 4278387459u, 4278453252u, 4278519045u, 4278584838u, 4278650631u, 4278716424u, 4278782217u, 4278848010u, 4278913803u, 4278979596u, 4279045389u, 4279111182u, 4279176975u, 4279242768u, 4279308561u, 4279374354u, 4279440147u, 4279505940u, 4279571733u, 4279637526u, 4279703319u, 4279769112u, 4279834905u, 4279900698u, 4279966491u, 4280032284u, 4280098077u, 4280163870u, 4280229663u, 4280295456u, 4280361249u, 4280427042u, 4280492835u, 4280558628u, 4280624421u, 4280690214u, 4280756007u, 4280821800u, 4280887593u, 4280953386u, 4281019179u, 4281084972u, 4281150765u, 4281216558u, 4281282351u, 4281348144u, 4281413937u, 4281479730u, 4281545523u, 4281611316u, 4281677109u, 4281742902u, 4281808695u, 4281874488u, 4281940281u, 4282006074u, 4282071867u, 4282137660u, 4282203453u, 4282269246u, 4282335039u, 4282400832u, 4282466625u, 4282532418u, 4282598211u, 4282664004u, 4282729797u, 4282795590u, 4282861383u, 4282927176u, 4282992969u, 4283058762u, 4283124555u, 4283190348u, 4283256141u, 4283321934u, 4283387727u, 4283453520u, 4283519313u, 4283585106u, 4283650899u, 4283716692u, 4283782485u, 4283848278u, 4283914071u, 4283979864u, 4284045657u, 4284111450u, 4284177243u, 4284243036u, 4284308829u, 4284374622u, 4284440415u, 4284506208u, 4284572001u, 4284637794u, 4284703587u, 4284769380u, 4284835173u, 4284900966u, 4284966759u, 4285032552u, 4285098345u, 4285164138u, 4285229931u, 4285295724u, 4285361517u, 4285427310u, 4285493103u, 4285558896u, 4285624689u, 4285690482u, 4285756275u, 4285822068u, 4285887861u, 4285953654u, 4286019447u, 4286085240u, 4286151033u, 4286216826u, 4286282619u, 4286348412u, 4286414205u, 4286479998u, 4286545791u, 4286611584u, 4286677377u, 4286743170u, 4286808963u, 4286874756u, 4286940549u, 4287006342u, 4287072135u, 4287137928u, 4287203721u, 4287269514u, 4287335307u, 4287401100u, 4287466893u, 4287532686u, 4287598479u, 4287664272u, 4287730065u, 4287795858u, 4287861651u, 4287927444u, 4287993237u, 4288059030u, 4288124823u, 4288190616u, 4288256409u, 4288322202u, 4288387995u, 4288453788u, 4288519581u, 4288585374u, 4288651167u, 4288716960u, 4288782753u, 4288848546u, 4288914339u, 4288980132u, 4289045925u, 4289111718u, 4289177511u, 4289243304u, 4289309097u, 4289374890u, 4289440683u, 4289506476u, 4289572269u, 4289638062u, 4289703855u, 4289769648u, 4289835441u, 4289901234u, 4289967027u, 4290032820u, 4290098613u, 4290164406u, 4290230199u, 4290295992u, 4290361785u, 4290427578u, 4290493371u, 4290559164u, 4290624957u, 4290690750u, 4290756543u, 4290822336u, 4290888129u, 4290953922u, 4291019715u, 4291085508u, 4291151301u, 4291217094u, 4291282887u, 4291348680u, 4291414473u, 4291480266u, 4291546059u, 4291611852u, 4291677645u, 4291743438u, 4291809231u, 4291875024u, 4291940817u, 4292006610u, 4292072403u, 4292138196u, 4292203989u, 4292269782u, 4292335575u, 4292401368u, 4292467161u, 4292532954u, 4292598747u, 4292664540u, 4292730333u, 4292796126u, 4292861919u, 4292927712u, 4292993505u, 4293059298u, 4293125091u, 4293190884u, 4293256677u, 4293322470u, 4293388263u, 4293454056u, 4293519849u, 4293585642u, 4293651435u, 4293717228u, 4293783021u, 4293848814u, 4293914607u, 4293980400u, 4294046193u, 4294111986u, 4294177779u, 4294243572u, 4294309365u, 4294375158u, 4294440951u, 4294506744u, 4294572537u, 4294638330u, 4294704123u, 4294769916u, 4294835709u, 4294901502u, 4294967295u
};

unsigned int __red2blue[256] = {
4294901760u, 4294902274u, 4294902788u, 4294903302u, 4294903816u, 4294904330u, 4294904844u, 4294905358u, 4294905872u, 4294906386u, 4294906900u, 4294907414u, 4294907928u, 4294908442u, 4294908956u, 4294909470u, 4294909984u, 4294910498u, 4294911012u, 4294911526u, 4294912040u, 4294912554u, 4294913068u, 4294913582u, 4294914096u, 4294914610u, 4294915124u, 4294915638u, 4294916152u, 4294916666u, 4294917180u, 4294917694u, 4294918208u, 4294918722u, 4294919236u, 4294919750u, 4294920264u, 4294920778u, 4294921292u, 4294921806u, 4294922320u, 4294922834u, 4294923348u, 4294923862u, 4294924376u, 4294924890u, 4294925404u, 4294925918u, 4294926432u, 4294926946u, 4294927460u, 4294927974u, 4294928488u, 4294929002u, 4294929516u, 4294930030u, 4294930544u, 4294931058u, 4294931572u, 4294932086u, 4294932600u, 4294933114u, 4294933628u, 4294934142u, 4294934656u, 4294935170u, 4294935684u, 4294936198u, 4294936712u, 4294937226u, 4294937740u, 4294938254u, 4294938768u, 4294939282u, 4294939796u, 4294940310u, 4294940824u, 4294941338u, 4294941852u, 4294942366u, 4294942880u, 4294943394u, 4294943908u, 4294944422u, 4294944936u, 4294945450u, 4294945964u, 4294946478u, 4294946992u, 4294947506u, 4294948020u, 4294948534u, 4294949048u, 4294949562u, 4294950076u, 4294950590u, 4294951104u, 4294951618u, 4294952132u, 4294952646u, 4294953160u, 4294953674u, 4294954188u, 4294954702u, 4294955216u, 4294955730u, 4294956244u, 4294956758u, 4294957272u, 4294957786u, 4294958300u, 4294958814u, 4294959328u, 4294959842u, 4294960356u, 4294960870u, 4294961384u, 4294961898u, 4294962412u, 4294962926u, 4294963440u, 4294963954u, 4294964468u, 4294964982u, 4294965496u, 4294966010u, 4294966524u, 4294967038u, 4294967295u, 4294835711u, 4294704127u, 4294572543u, 4294440959u, 4294309375u, 4294177791u, 4294046207u, 4293914623u, 4293783039u, 4293651455u, 4293519871u, 4293388287u, 4293256703u, 4293125119u, 4292993535u, 4292861951u, 4292730367u, 4292598783u, 4292467199u, 4292335615u, 4292204031u, 4292072447u, 4291940863u, 4291809279u, 4291677695u, 4291546111u, 4291414527u, 4291282943u, 4291151359u, 4291019775u, 4290888191u, 4290756607u, 4290625023u, 4290493439u, 4290361855u, 4290230271u, 4290098687u, 4289967103u, 4289835519u, 4289703935u, 4289572351u, 4289440767u, 4289309183u, 4289177599u, 4289046015u, 4288914431u, 4288782847u, 4288651263u, 4288519679u, 4288388095u, 4288256511u, 4288124927u, 4287993343u, 4287861759u, 4287730175u, 4287598591u, 4287467007u, 4287335423u, 4287203839u, 4287072255u, 4286940671u, 4286809087u, 4286677503u, 4286545919u, 4286414335u, 4286282751u, 4286151167u, 4286019583u, 4285887999u, 4285756415u, 4285624831u, 4285493247u, 4285361663u, 4285230079u, 4285098495u, 4284966911u, 4284835327u, 4284703743u, 4284572159u, 4284440575u, 4284308991u, 4284177407u, 4284045823u, 4283914239u, 4283782655u, 4283651071u, 4283519487u, 4283387903u, 4283256319u, 4283124735u, 4282993151u, 4282861567u, 4282729983u, 4282598399u, 4282466815u, 4282335231u, 4282203647u, 4282072063u, 4281940479u, 4281808895u, 4281677311u, 4281545727u, 4281414143u, 4281282559u, 4281150975u, 4281019391u, 4280887807u, 4280756223u, 4280624639u, 4280493055u, 4280361471u, 4280229887u, 4280098303u, 4279966719u, 4279835135u, 4279703551u, 4279571967u, 4279440383u, 4279308799u, 4279177215u, 4279045631u, 4278914047u, 4278782463u, 4278650879u, 4278519295u, 4278387711u, 4278256127u
};

unsigned int __blue2red[256] = {
    4278190335u, 4278321919u, 4278453503u, 4278585087u, 4278716671u, 4278848255u, 4278979839u, 4279111423u, 4279243007u, 4279374591u, 4279506175u, 4279637759u, 4279769343u, 4279900927u, 4280032511u, 4280164095u, 4280295679u, 4280427263u, 4280558847u, 4280690431u, 4280822015u, 4280953599u, 4281085183u, 4281216767u, 4281348351u, 4281479935u, 4281611519u, 4281743103u, 4281874687u, 4282006271u, 4282137855u, 4282269439u, 4282401023u, 4282532607u, 4282664191u, 4282795775u, 4282927359u, 4283058943u, 4283190527u, 4283322111u, 4283453695u, 4283585279u, 4283716863u, 4283848447u, 4283980031u, 4284111615u, 4284243199u, 4284374783u, 4284506367u, 4284637951u, 4284769535u, 4284901119u, 4285032703u, 4285164287u, 4285295871u, 4285427455u, 4285559039u, 4285690623u, 4285822207u, 4285953791u, 4286085375u, 4286216959u, 4286348543u, 4286480127u, 4286611711u, 4286743295u, 4286874879u, 4287006463u, 4287138047u, 4287269631u, 4287401215u, 4287532799u, 4287664383u, 4287795967u, 4287927551u, 4288059135u, 4288190719u, 4288322303u, 4288453887u, 4288585471u, 4288717055u, 4288848639u, 4288980223u, 4289111807u, 4289243391u, 4289374975u, 4289506559u, 4289638143u, 4289769727u, 4289901311u, 4290032895u, 4290164479u, 4290296063u, 4290427647u, 4290559231u, 4290690815u, 4290822399u, 4290953983u, 4291085567u, 4291217151u, 4291348735u, 4291480319u, 4291611903u, 4291743487u, 4291875071u, 4292006655u, 4292138239u, 4292269823u, 4292401407u, 4292532991u, 4292664575u, 4292796159u, 4292927743u, 4293059327u, 4293190911u, 4293322495u, 4293454079u, 4293585663u, 4293717247u, 4293848831u, 4293980415u, 4294111999u, 4294243583u, 4294375167u, 4294506751u, 4294638335u, 4294769919u, 4294901503u, 4294967295u, 4294966781u, 4294966267u, 4294965753u, 4294965239u, 4294964725u, 4294964211u, 4294963697u, 4294963183u, 4294962669u, 4294962155u, 4294961641u, 4294961127u, 4294960613u, 4294960099u, 4294959585u, 4294959071u, 4294958557u, 4294958043u, 4294957529u, 4294957015u, 4294956501u, 4294955987u, 4294955473u, 4294954959u, 4294954445u, 4294953931u, 4294953417u, 4294952903u, 4294952389u, 4294951875u, 4294951361u, 4294950847u, 4294950333u, 4294949819u, 4294949305u, 4294948791u, 4294948277u, 4294947763u, 4294947249u, 4294946735u, 4294946221u, 4294945707u, 4294945193u, 4294944679u, 4294944165u, 4294943651u, 4294943137u, 4294942623u, 4294942109u, 4294941595u, 4294941081u, 4294940567u, 4294940053u, 4294939539u, 4294939025u, 4294938511u, 4294937997u, 4294937483u, 4294936969u, 4294936455u, 4294935941u, 4294935427u, 4294934913u, 4294934399u, 4294933885u, 4294933371u, 4294932857u, 4294932343u, 4294931829u, 4294931315u, 4294930801u, 4294930287u, 4294929773u, 4294929259u, 4294928745u, 4294928231u, 4294927717u, 4294927203u, 4294926689u, 4294926175u, 4294925661u, 4294925147u, 4294924633u, 4294924119u, 4294923605u, 4294923091u, 4294922577u, 4294922063u, 4294921549u, 4294921035u, 4294920521u, 4294920007u, 4294919493u, 4294918979u, 4294918465u, 4294917951u, 4294917437u, 4294916923u, 4294916409u, 4294915895u, 4294915381u, 4294914867u, 4294914353u, 4294913839u, 4294913325u, 4294912811u, 4294912297u, 4294911783u, 4294911269u, 4294910755u, 4294910241u, 4294909727u, 4294909213u, 4294908699u, 4294908185u, 4294907671u, 4294907157u, 4294906643u, 4294906129u, 4294905615u, 4294905101u, 4294904587u, 4294904073u, 4294903559u, 4294903045u, 4294902531u, 4294902017u
};

unsigned int __blue2redX[256] = {
    4278190335u, 4278321919u, 4278453503u, 4278585087u, 4278716671u, 4278848255u, 4278979839u, 4279111423u, 4279243007u, 4279374591u, 4279506175u, 4279637759u, 4279769343u, 4279900927u, 4280032511u, 4280164095u, 4280295679u, 4280427263u, 4280558847u, 4280690431u, 4280822015u, 4280953599u, 4281085183u, 4281216767u, 4281348351u, 4281479935u, 4281611519u, 4281743103u, 4281874687u, 4282006271u, 4282137855u, 4282269439u, 4282401023u, 4282532607u, 4282664191u, 4282795775u, 4282927359u, 4283058943u, 4283190527u, 4283322111u, 4283453695u, 4283585279u, 4283716863u, 4283848447u, 4283980031u, 4284111615u, 4284243199u, 4284374783u, 4284506367u, 4284637951u, 4284769535u, 4284901119u, 4285032703u, 4285164287u, 4285295871u, 4285427455u, 4285559039u, 4285690623u, 4285822207u, 4285953791u, 4286085375u, 4286216959u, 4286348543u, 4286480127u, 4286611711u, 4286743295u, 4286874879u, 4287006463u, 4287138047u, 4287269631u, 4287401215u, 4287532799u, 4287664383u, 4287795967u, 4287927551u, 4288059135u, 4288190719u, 4288322303u, 4288453887u, 4288585471u, 4288717055u, 4288848639u, 4288980223u, 4289111807u, 4289243391u, 4289374975u, 4289506559u, 4289638143u, 4289769727u, 4289901311u, 4290032895u, 4290164479u, 4290296063u, 4290427647u, 4290559231u, 4290690815u, 4290822399u, 4290953983u, 4291085567u, 4291217151u, 4291348735u, 4291480319u, 4291611903u, 4291743487u, 4291875071u, 4292006655u, 4292138239u, 4292269823u, 4292401407u, 4292532991u, 4292664575u, 4292796159u, 4292927743u, 4293059327u, 4293190911u, 4293322495u, 4293454079u, 4293585663u, 4293717247u, 4293848831u, 4293980415u, 4294111999u, 4294243583u, 4294375167u, 4294506751u, 4294638335u, 4294769919u, 4294901503u, 4294967295u, 4294966781u, 4294966267u, 4294965753u, 4294965239u, 4294964725u, 4294964211u, 4294963697u, 4294963183u, 4294962669u, 4294962155u, 4294961641u, 4294961127u, 4294960613u, 4294960099u, 4294959585u, 4294959071u, 4294958557u, 4294958043u, 4294957529u, 4294957015u, 4294956501u, 4294955987u, 4294955473u, 4294954959u, 4294954445u, 4294953931u, 4294953417u, 4294952903u, 4294952389u, 4294951875u, 4294951361u, 4294950847u, 4294950333u, 4294949819u, 4294949305u, 4294948791u, 4294948277u, 4294947763u, 4294947249u, 4294946735u, 4294946221u, 4294945707u, 4294945193u, 4294944679u, 4294944165u, 4294943651u, 4294943137u, 4294942623u, 4294942109u, 4294941595u, 4294941081u, 4294940567u, 4294940053u, 4294939539u, 4294939025u, 4294938511u, 4294937997u, 4294937483u, 4294936969u, 4294936455u, 4294935941u, 4294935427u, 4294934913u, 4294934399u, 4294933885u, 4294933371u, 4294932857u, 4294932343u, 4294931829u, 4294931315u, 4294930801u, 4294930287u, 4294929773u, 4294929259u, 4294928745u, 4294928231u, 4294927717u, 4294927203u, 4294926689u, 4294926175u, 4294925661u, 4294925147u, 4294924633u, 4294924119u, 4294923605u, 4294923091u, 4294922577u, 4294922063u, 4294921549u, 4294921035u, 4294920521u, 4294920007u, 4294919493u, 4294918979u, 4294918465u, 4294917951u, 4294917437u, 4294916923u, 4294916409u, 4294915895u, 4294915381u, 4294914867u, 4294914353u, 4294913839u, 4294913325u, 4294912811u, 4294912297u, 4294911783u, 4294911269u, 4294910755u, 4294910241u, 4294909727u, 4294909213u, 4294908699u, 4294908185u, 4294907671u, 4294907157u, 4294906643u, 4294906129u, 4294905615u, 4294905101u, 4294904587u, 4294904073u, 4294903559u, 4294903045u, 4294902531u, 0
};

unsigned int __hsv1000[1000] = {
    0xffff0000,0xffff0100,0xffff0300,0xffff0400,0xffff0600,0xffff0700,0xffff0900,0xffff0a00,0xffff0c00,0xffff0d00,0xffff0f00,0xffff1000,0xffff1200,0xffff1300,0xffff1500,0xffff1600,0xffff1800,0xffff1a00,0xffff1b00,0xffff1d00,0xffff1e00,0xffff2000,0xffff2100,0xffff2300,0xffff2400,0xffff2600,0xffff2700,0xffff2900,0xffff2a00,0xffff2c00,0xffff2d00,0xffff2f00,0xffff3000,0xffff3200,0xffff3400,0xffff3500,0xffff3700,0xffff3800,0xffff3a00,0xffff3b00,0xffff3d00,0xffff3e00,0xffff4000,0xffff4100,0xffff4300,0xffff4400,0xffff4600,0xffff4700,0xffff4900,0xffff4a00,0xffff4c00,0xffff4e00,0xffff4f00,0xffff5100,0xffff5200,0xffff5400,0xffff5500,0xffff5700,0xffff5800,0xffff5a00,0xffff5b00,0xffff5d00,0xffff5e00,0xffff6000,0xffff6100,0xffff6300,0xffff6400,0xffff6600,0xffff6800,0xffff6900,0xffff6b00,0xffff6c00,0xffff6e00,0xffff6f00,0xffff7100,0xffff7200,0xffff7400,0xffff7500,0xffff7700,0xffff7800,0xffff7a00,0xffff7b00,0xffff7d00,0xffff7e00,0xffff8000,0xffff8200,0xffff8300,0xffff8500,0xffff8600,0xffff8800,0xffff8900,0xffff8b00,0xffff8c00,0xffff8e00,0xffff8f00,0xffff9100,0xffff9200,0xffff9400,0xffff9500,0xffff9700,0xffff9900,0xffff9a00,0xffff9c00,0xffff9d00,0xffff9f00,0xffffa000,0xffffa200,0xffffa300,0xffffa500,0xffffa600,0xffffa800,0xffffa900,0xffffab00,0xffffac00,0xffffae00,0xffffaf00,0xffffb100,0xffffb300,0xffffb400,0xffffb600,0xffffb700,0xffffb900,0xffffba00,0xffffbc00,0xffffbd00,0xffffbf00,0xffffc000,0xffffc200,0xffffc300,0xffffc500,0xffffc600,0xffffc800,0xffffc900,0xffffcb00,0xffffcd00,0xffffce00,0xffffd000,0xffffd100,0xffffd300,0xffffd400,0xffffd600,0xffffd700,0xffffd900,0xffffda00,0xffffdc00,0xffffdd00,0xffffdf00,0xffffe000,0xffffe200,0xffffe300,0xffffe500,0xffffe700,0xffffe800,0xffffea00,0xffffeb00,0xffffed00,0xffffee00,0xfffff000,0xfffff100,0xfffff300,0xfffff400,0xfffff600,0xfffff700,0xfffff900,0xfffffa00,0xfffffc00,0xfffffd00,0xfffeff00,0xfffcff00,0xfffbff00,0xfff9ff00,0xfff8ff00,0xfff6ff00,0xfff5ff00,0xfff3ff00,0xfff2ff00,0xfff0ff00,0xffefff00,0xffedff00,0xffecff00,0xffeaff00,0xffe9ff00,0xffe7ff00,0xffe6ff00,0xffe4ff00,0xffe2ff00,0xffe1ff00,0xffdfff00,0xffdeff00,0xffdcff00,0xffdbff00,0xffd9ff00,0xffd8ff00,0xffd6ff00,0xffd5ff00,0xffd3ff00,0xffd2ff00,0xffd0ff00,0xffcfff00,0xffcdff00,0xffcbff00,0xffcaff00,0xffc8ff00,0xffc7ff00,0xffc5ff00,0xffc4ff00,0xffc2ff00,0xffc1ff00,0xffbfff00,0xffbeff00,0xffbcff00,0xffbbff00,0xffb9ff00,0xffb8ff00,0xffb6ff00,0xffb5ff00,0xffb3ff00,0xffb1ff00,0xffb0ff00,0xffaeff00,0xffadff00,0xffabff00,0xffaaff00,0xffa8ff00,0xffa7ff00,0xffa5ff00,0xffa4ff00,0xffa2ff00,0xffa1ff00,0xff9fff00,0xff9eff00,0xff9cff00,0xff9bff00,0xff99ff00,0xff97ff00,0xff96ff00,0xff94ff00,0xff93ff00,0xff91ff00,0xff90ff00,0xff8eff00,0xff8dff00,0xff8bff00,0xff8aff00,0xff88ff00,0xff87ff00,0xff85ff00,0xff84ff00,0xff82ff00,0xff81ff00,0xff7fff00,0xff7dff00,0xff7cff00,0xff7aff00,0xff79ff00,0xff77ff00,0xff76ff00,0xff74ff00,0xff73ff00,0xff71ff00,0xff70ff00,0xff6eff00,0xff6dff00,0xff6bff00,0xff6aff00,0xff68ff00,0xff67ff00,0xff65ff00,0xff63ff00,0xff62ff00,0xff60ff00,0xff5fff00,0xff5dff00,0xff5cff00,0xff5aff00,0xff59ff00,0xff57ff00,0xff56ff00,0xff54ff00,0xff53ff00,0xff51ff00,0xff50ff00,0xff4eff00,0xff4dff00,0xff4bff00,0xff49ff00,0xff48ff00,0xff46ff00,0xff45ff00,0xff43ff00,0xff42ff00,0xff40ff00,0xff3fff00,0xff3dff00,0xff3cff00,0xff3aff00,0xff39ff00,0xff37ff00,0xff36ff00,0xff34ff00,0xff33ff00,0xff31ff00,0xff2fff00,0xff2eff00,0xff2cff00,0xff2bff00,0xff29ff00,0xff28ff00,0xff26ff00,0xff25ff00,0xff23ff00,0xff22ff00,0xff20ff00,0xff1fff00,0xff1dff00,0xff1cff00,0xff1aff00,0xff18ff00,0xff17ff00,0xff15ff00,0xff14ff00,0xff12ff00,0xff11ff00,0xff0fff00,0xff0eff00,0xff0cff00,0xff0bff00,0xff09ff00,0xff08ff00,0xff06ff00,0xff05ff00,0xff03ff00,0xff02ff00,0xff00ff00,0xff00ff01,0xff00ff02,0xff00ff04,0xff00ff05,0xff00ff07,0xff00ff08,0xff00ff0a,0xff00ff0b,0xff00ff0d,0xff00ff0e,0xff00ff10,0xff00ff11,0xff00ff13,0xff00ff14,0xff00ff16,0xff00ff17,0xff00ff19,0xff00ff1b,0xff00ff1c,0xff00ff1e,0xff00ff1f,0xff00ff21,0xff00ff22,0xff00ff24,0xff00ff25,0xff00ff27,0xff00ff28,0xff00ff2a,0xff00ff2b,0xff00ff2d,0xff00ff2e,0xff00ff30,0xff00ff31,0xff00ff33,0xff00ff35,0xff00ff36,0xff00ff38,0xff00ff39,0xff00ff3b,0xff00ff3c,0xff00ff3e,0xff00ff3f,0xff00ff41,0xff00ff42,0xff00ff44,0xff00ff45,0xff00ff47,0xff00ff48,0xff00ff4a,0xff00ff4b,0xff00ff4d,0xff00ff4f,0xff00ff50,0xff00ff52,0xff00ff53,0xff00ff55,0xff00ff56,0xff00ff58,0xff00ff59,0xff00ff5b,0xff00ff5c,0xff00ff5e,0xff00ff5f,0xff00ff61,0xff00ff62,0xff00ff64,0xff00ff66,0xff00ff67,0xff00ff69,0xff00ff6a,0xff00ff6c,0xff00ff6d,0xff00ff6f,0xff00ff70,0xff00ff72,0xff00ff73,0xff00ff75,0xff00ff76,0xff00ff78,0xff00ff79,0xff00ff7b,0xff00ff7c,0xff00ff7e,0xff00ff80,0xff00ff81,0xff00ff83,0xff00ff84,0xff00ff86,0xff00ff87,0xff00ff89,0xff00ff8a,0xff00ff8c,0xff00ff8d,0xff00ff8f,0xff00ff90,0xff00ff92,0xff00ff93,0xff00ff95,0xff00ff96,0xff00ff98,0xff00ff9a,0xff00ff9b,0xff00ff9d,0xff00ff9e,0xff00ffa0,0xff00ffa1,0xff00ffa3,0xff00ffa4,0xff00ffa6,0xff00ffa7,0xff00ffa9,0xff00ffaa,0xff00ffac,0xff00ffad,0xff00ffaf,0xff00ffb0,0xff00ffb2,0xff00ffb4,0xff00ffb5,0xff00ffb7,0xff00ffb8,0xff00ffba,0xff00ffbb,0xff00ffbd,0xff00ffbe,0xff00ffc0,0xff00ffc1,0xff00ffc3,0xff00ffc4,0xff00ffc6,0xff00ffc7,0xff00ffc9,0xff00ffca,0xff00ffcc,0xff00ffce,0xff00ffcf,0xff00ffd1,0xff00ffd2,0xff00ffd4,0xff00ffd5,0xff00ffd7,0xff00ffd8,0xff00ffda,0xff00ffdb,0xff00ffdd,0xff00ffde,0xff00ffe0,0xff00ffe1,0xff00ffe3,0xff00ffe4,0xff00ffe6,0xff00ffe8,0xff00ffe9,0xff00ffeb,0xff00ffec,0xff00ffee,0xff00ffef,0xff00fff1,0xff00fff2,0xff00fff4,0xff00fff5,0xff00fff7,0xff00fff8,0xff00fffa,0xff00fffb,0xff00fffd,0xff00ffff,0xff00fdff,0xff00fbff,0xff00faff,0xff00f8ff,0xff00f7ff,0xff00f5ff,0xff00f4ff,0xff00f2ff,0xff00f1ff,0xff00efff,0xff00eeff,0xff00ecff,0xff00ebff,0xff00e9ff,0xff00e8ff,0xff00e6ff,0xff00e4ff,0xff00e3ff,0xff00e1ff,0xff00e0ff,0xff00deff,0xff00ddff,0xff00dbff,0xff00daff,0xff00d8ff,0xff00d7ff,0xff00d5ff,0xff00d4ff,0xff00d2ff,0xff00d1ff,0xff00cfff,0xff00ceff,0xff00ccff,0xff00caff,0xff00c9ff,0xff00c7ff,0xff00c6ff,0xff00c4ff,0xff00c3ff,0xff00c1ff,0xff00c0ff,0xff00beff,0xff00bdff,0xff00bbff,0xff00baff,0xff00b8ff,0xff00b7ff,0xff00b5ff,0xff00b4ff,0xff00b2ff,0xff00b0ff,0xff00afff,0xff00adff,0xff00acff,0xff00aaff,0xff00a9ff,0xff00a7ff,0xff00a6ff,0xff00a4ff,0xff00a3ff,0xff00a1ff,0xff00a0ff,0xff009eff,0xff009dff,0xff009bff,0xff009aff,0xff0098ff,0xff0096ff,0xff0095ff,0xff0093ff,0xff0092ff,0xff0090ff,0xff008fff,0xff008dff,0xff008cff,0xff008aff,0xff0089ff,0xff0087ff,0xff0086ff,0xff0084ff,0xff0083ff,0xff0081ff,0xff0080ff,0xff007eff,0xff007cff,0xff007bff,0xff0079ff,0xff0078ff,0xff0076ff,0xff0075ff,0xff0073ff,0xff0072ff,0xff0070ff,0xff006fff,0xff006dff,0xff006cff,0xff006aff,0xff0069ff,0xff0067ff,0xff0066ff,0xff0064ff,0xff0062ff,0xff0061ff,0xff005fff,0xff005eff,0xff005cff,0xff005bff,0xff0059ff,0xff0058ff,0xff0056ff,0xff0055ff,0xff0053ff,0xff0052ff,0xff0050ff,0xff004fff,0xff004dff,0xff004bff,0xff004aff,0xff0048ff,0xff0047ff,0xff0045ff,0xff0044ff,0xff0042ff,0xff0041ff,0xff003fff,0xff003eff,0xff003cff,0xff003bff,0xff0039ff,0xff0038ff,0xff0036ff,0xff0035ff,0xff0033ff,0xff0031ff,0xff0030ff,0xff002eff,0xff002dff,0xff002bff,0xff002aff,0xff0028ff,0xff0027ff,0xff0025ff,0xff0024ff,0xff0022ff,0xff0021ff,0xff001fff,0xff001eff,0xff001cff,0xff001bff,0xff0019ff,0xff0017ff,0xff0016ff,0xff0014ff,0xff0013ff,0xff0011ff,0xff0010ff,0xff000eff,0xff000dff,0xff000bff,0xff000aff,0xff0008ff,0xff0007ff,0xff0005ff,0xff0004ff,0xff0002ff,0xff0001ff,0xff0000ff,0xff0200ff,0xff0300ff,0xff0500ff,0xff0600ff,0xff0800ff,0xff0900ff,0xff0b00ff,0xff0c00ff,0xff0e00ff,0xff0f00ff,0xff1100ff,0xff1200ff,0xff1400ff,0xff1500ff,0xff1700ff,0xff1800ff,0xff1a00ff,0xff1c00ff,0xff1d00ff,0xff1f00ff,0xff2000ff,0xff2200ff,0xff2300ff,0xff2500ff,0xff2600ff,0xff2800ff,0xff2900ff,0xff2b00ff,0xff2c00ff,0xff2e00ff,0xff2f00ff,0xff3100ff,0xff3300ff,0xff3400ff,0xff3600ff,0xff3700ff,0xff3900ff,0xff3a00ff,0xff3c00ff,0xff3d00ff,0xff3f00ff,0xff4000ff,0xff4200ff,0xff4300ff,0xff4500ff,0xff4600ff,0xff4800ff,0xff4900ff,0xff4b00ff,0xff4d00ff,0xff4e00ff,0xff5000ff,0xff5100ff,0xff5300ff,0xff5400ff,0xff5600ff,0xff5700ff,0xff5900ff,0xff5a00ff,0xff5c00ff,0xff5d00ff,0xff5f00ff,0xff6000ff,0xff6200ff,0xff6300ff,0xff6500ff,0xff6700ff,0xff6800ff,0xff6a00ff,0xff6b00ff,0xff6d00ff,0xff6e00ff,0xff7000ff,0xff7100ff,0xff7300ff,0xff7400ff,0xff7600ff,0xff7700ff,0xff7900ff,0xff7a00ff,0xff7c00ff,0xff7d00ff,0xff7f00ff,0xff8100ff,0xff8200ff,0xff8400ff,0xff8500ff,0xff8700ff,0xff8800ff,0xff8a00ff,0xff8b00ff,0xff8d00ff,0xff8e00ff,0xff9000ff,0xff9100ff,0xff9300ff,0xff9400ff,0xff9600ff,0xff9700ff,0xff9900ff,0xff9b00ff,0xff9c00ff,0xff9e00ff,0xff9f00ff,0xffa100ff,0xffa200ff,0xffa400ff,0xffa500ff,0xffa700ff,0xffa800ff,0xffaa00ff,0xffab00ff,0xffad00ff,0xffae00ff,0xffb000ff,0xffb100ff,0xffb300ff,0xffb500ff,0xffb600ff,0xffb800ff,0xffb900ff,0xffbb00ff,0xffbc00ff,0xffbe00ff,0xffbf00ff,0xffc100ff,0xffc200ff,0xffc400ff,0xffc500ff,0xffc700ff,0xffc800ff,0xffca00ff,0xffcc00ff,0xffcd00ff,0xffcf00ff,0xffd000ff,0xffd200ff,0xffd300ff,0xffd500ff,0xffd600ff,0xffd800ff,0xffd900ff,0xffdb00ff,0xffdc00ff,0xffde00ff,0xffdf00ff,0xffe100ff,0xffe200ff,0xffe400ff,0xffe600ff,0xffe700ff,0xffe900ff,0xffea00ff,0xffec00ff,0xffed00ff,0xffef00ff,0xfff000ff,0xfff200ff,0xfff300ff,0xfff500ff,0xfff600ff,0xfff800ff,0xfff900ff,0xfffb00ff,0xfffc00ff,0xfffe00ff,0xffff00fd,0xffff00fc,0xffff00fa,0xffff00f9,0xffff00f7,0xffff00f6,0xffff00f4,0xffff00f3,0xffff00f1,0xffff00f0,0xffff00ee,0xffff00ed,0xffff00eb,0xffff00ea,0xffff00e8,0xffff00e7,0xffff00e5,0xffff00e3,0xffff00e2,0xffff00e0,0xffff00df,0xffff00dd,0xffff00dc,0xffff00da,0xffff00d9,0xffff00d7,0xffff00d6,0xffff00d4,0xffff00d3,0xffff00d1,0xffff00d0,0xffff00ce,0xffff00cd,0xffff00cb,0xffff00c9,0xffff00c8,0xffff00c6,0xffff00c5,0xffff00c3,0xffff00c2,0xffff00c0,0xffff00bf,0xffff00bd,0xffff00bc,0xffff00ba,0xffff00b9,0xffff00b7,0xffff00b6,0xffff00b4,0xffff00b3,0xffff00b1,0xffff00af,0xffff00ae,0xffff00ac,0xffff00ab,0xffff00a9,0xffff00a8,0xffff00a6,0xffff00a5,0xffff00a3,0xffff00a2,0xffff00a0,0xffff009f,0xffff009d,0xffff009c,0xffff009a,0xffff0098,0xffff0097,0xffff0095,0xffff0094,0xffff0092,0xffff0091,0xffff008f,0xffff008e,0xffff008c,0xffff008b,0xffff0089,0xffff0088,0xffff0086,0xffff0085,0xffff0083,0xffff0082,0xffff0080,0xffff007e,0xffff007d,0xffff007b,0xffff007a,0xffff0078,0xffff0077,0xffff0075,0xffff0074,0xffff0072,0xffff0071,0xffff006f,0xffff006e,0xffff006c,0xffff006b,0xffff0069,0xffff0068,0xffff0066,0xffff0064,0xffff0063,0xffff0061,0xffff0060,0xffff005e,0xffff005d,0xffff005b,0xffff005a,0xffff0058,0xffff0057,0xffff0055,0xffff0054,0xffff0052,0xffff0051,0xffff004f,0xffff004e,0xffff004c,0xffff004a,0xffff0049,0xffff0047,0xffff0046,0xffff0044,0xffff0043,0xffff0041,0xffff0040,0xffff003e,0xffff003d,0xffff003b,0xffff003a,0xffff0038,0xffff0037,0xffff0035,0xffff0034,0xffff0032,0xffff0030,0xffff002f,0xffff002d,0xffff002c,0xffff002a,0xffff0029,0xffff0027,0xffff0026,0xffff0024,0xffff0023,0xffff0021,0xffff0020,0xffff001e,0xffff001d,0xffff001b,0xffff001a,0xffff0018,0xffff0016,0xffff0015,0xffff0013,0xffff0012,0xffff0010,0xffff000f,0xffff000d,0xffff000c,0xffff000a,0xffff0009,0xffff0007,0xffff0006,0xffff0004,0xffff0003,0xff5555555//,0xffff0001
};

unsigned int __hsv256[256] = {
    0xffff0000,0xffff0005,0xffff000b,0xffff0011,0xffff0017,0xffff001d,0xffff0023,0xffff0029,0xffff002f,0xffff0035,0xffff003b,0xffff0041,0xffff0047,0xffff004d,0xffff0053,0xffff0059,0xffff005f,0xffff0065,0xffff006b,0xffff0071,0xffff0077,0xffff007d,0xffff0083,0xffff0089,0xffff008f,0xffff0095,0xffff009b,0xffff00a1,0xffff00a7,0xffff00ad,0xffff00b3,0xffff00b9,0xffff00bf,0xffff00c5,0xffff00cb,0xffff00d1,0xffff00d7,0xffff00dd,0xffff00e3,0xffff00e9,0xffff00ef,0xffff00f5,0xffff00fb,0xfffd00ff,0xfff700ff,0xfff100ff,0xffeb00ff,0xffe500ff,0xffdf00ff,0xffd900ff,0xffd300ff,0xffcd00ff,0xffc700ff,0xffc100ff,0xffbb00ff,0xffb500ff,0xffaf00ff,0xffa900ff,0xffa300ff,0xff9d00ff,0xff9700ff,0xff9100ff,0xff8b00ff,0xff8500ff,0xff7f00ff,0xff7900ff,0xff7300ff,0xff6d00ff,0xff6700ff,0xff6100ff,0xff5b00ff,0xff5500ff,0xff4f00ff,0xff4900ff,0xff4300ff,0xff3d00ff,0xff3700ff,0xff3100ff,0xff2b00ff,0xff2500ff,0xff1f00ff,0xff1900ff,0xff1300ff,0xff0d00ff,0xff0700ff,0xff0100ff,0xff0003ff,0xff0009ff,0xff000fff,0xff0015ff,0xff001bff,0xff0021ff,0xff0027ff,0xff002dff,0xff0033ff,0xff0039ff,0xff003fff,0xff0045ff,0xff004bff,0xff0051ff,0xff0057ff,0xff005dff,0xff0063ff,0xff0069ff,0xff006fff,0xff0075ff,0xff007bff,0xff0081ff,0xff0087ff,0xff008dff,0xff0093ff,0xff0099ff,0xff009fff,0xff00a5ff,0xff00abff,0xff00b1ff,0xff00b7ff,0xff00bdff,0xff00c3ff,0xff00c9ff,0xff00cfff,0xff00d5ff,0xff00dbff,0xff00e1ff,0xff00e7ff,0xff00edff,0xff00f3ff,0xff00f9ff,0xff00ffff,0xff00fff9,0xff00fff3,0xff00ffed,0xff00ffe7,0xff00ffe1,0xff00ffdb,0xff00ffd5,0xff00ffcf,0xff00ffc9,0xff00ffc3,0xff00ffbd,0xff00ffb7,0xff00ffb1,0xff00ffab,0xff00ffa5,0xff00ff9f,0xff00ff99,0xff00ff93,0xff00ff8d,0xff00ff87,0xff00ff81,0xff00ff7b,0xff00ff75,0xff00ff6f,0xff00ff69,0xff00ff63,0xff00ff5d,0xff00ff57,0xff00ff51,0xff00ff4b,0xff00ff45,0xff00ff3f,0xff00ff39,0xff00ff33,0xff00ff2d,0xff00ff27,0xff00ff21,0xff00ff1b,0xff00ff15,0xff00ff0f,0xff00ff09,0xff00ff03,0xff01ff00,0xff07ff00,0xff0dff00,0xff13ff00,0xff19ff00,0xff1fff00,0xff25ff00,0xff2bff00,0xff31ff00,0xff37ff00,0xff3dff00,0xff43ff00,0xff49ff00,0xff4fff00,0xff55ff00,0xff5bff00,0xff61ff00,0xff67ff00,0xff6dff00,0xff73ff00,0xff79ff00,0xff7fff00,0xff85ff00,0xff8bff00,0xff91ff00,0xff97ff00,0xff9dff00,0xffa3ff00,0xffa9ff00,0xffafff00,0xffb5ff00,0xffbbff00,0xffc1ff00,0xffc7ff00,0xffcdff00,0xffd3ff00,0xffd9ff00,0xffdfff00,0xffe5ff00,0xffebff00,0xfff1ff00,0xfff7ff00,0xfffdff00,0xfffffb00,0xfffff500,0xffffef00,0xffffe900,0xffffe300,0xffffdd00,0xffffd700,0xffffd100,0xffffcb00,0xffffc500,0xffffbf00,0xffffb900,0xffffb300,0xffffad00,0xffffa700,0xffffa100,0xffff9b00,0xffff9500,0xffff8f00,0xffff8900,0xffff8300,0xffff7d00,0xffff7700,0xffff7100,0xffff6b00,0xffff6500,0xffff5f00,0xffff5900,0xffff5300,0xffff4d00,0xffff4700,0xffff4100,0xffff3b00,0xffff3500,0xffff2f00,0xffff2900,0xffff2300,0xffff1d00,0xffff1700,0xffff1100,0xffff0b00,0xff666666//0xffff0500
};