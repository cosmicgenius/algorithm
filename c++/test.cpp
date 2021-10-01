#include "lib/Algorithm.h"
#include "lib/Bint.h"

using namespace std;

int main() {
  // freopen("test.txt", "w", stdout);

  cout << "start" << endl;
  // Bint b1 =
  //     "116650792548111340344612018800484293276843391972748158301915781972989260"
  //     "432243812959272959023682044933847714492379681406129049042740760037008001"
  //     "354678327854523986007875892024005922406944891211108898373793144414776514"
  //     "512660601042546145692359791982382978018729020879606557621807853093592682"
  //     "868177207606816098751632428927777360329407982961867671548991418062873170"
  //     "829960930139503752725722683502664009625451186661321162099457536657840003"
  //     "123629805598332698460649794269581920378336688302063725165260485075892826"
  //     "464857860292283300233207249511405453494520805850905796263002377749078829"
  //     "799841583043754618844158819675002093888832325483313596179094599653245083"
  //     "810270180137214271274145880682366495331083889820174101073843265838118853"
  //     "779574725310147436979636025215268722598714017667995069286215325039584252"
  //     "781545961725334010503932750350772209940788987590223088368806937575384551"
  //     "319329381251368406532213398840980361905430203423522810322957938432099644"
  //     "056523994680538670255092095247382681253596261352370863840571082990872889"
  //     "171167861808569166312564616097314195998232073432863782510304892196896481"
  //     "955953213817453037568623165547567755102648889260800495578108283744048148"
  //     "770601497889367482665328551990356779453880357961846778577336044646675630"
  //     "588458772280224922417436340479660165960734464537627898443726354639290104"
  //     "972866941136790019499234865592487615032184517021729314476086666733714822"
  //     "501969466451513725035330353366403";
  // Bint b2 =
  //     "795638739407130611865013118588335678167919732446402422306806800186333102"
  //     "468966084096066976429285748690498343019527733504019035541614547845481862"
  //     "910130953201954574460432020976672571796406250443987900743008586349701194"
  //     "342767927822311057469577054452996489360825023025860802060619530963692827"
  //     "873161475519822390431048293639742766351747872271344783642561547539790338"
  //     "837273770106380956580663196769863922619314510491377116536200414343659777"
  //     "686859763378112290820760579146012338547407069570772998675070756872380635"
  //     "380762619622858005901191125324953029211301339300325034604657569486151870"
  //     "055082316700844659487081174715767851725187801683866132684988941873778948"
  //     "084827867508944348070072338467373420949683680201360970034868714699081749"
  //     "206150451536005333081536345475981293835362824834932799994341948679977874"
  //     "074331468590507837981593106725920956725475184060230170442048255670840115"
  //     "508877029107493692395905980886587798473240373718391362901549342955475521"
  //     "293890373952667927746788333108266308001905543326794665523704120733380228"
  //     "608370932627153006544178117376750944485690418437888355326896088818821584"
  //     "923091607940466120963119109397031225833388391618490992343753747618114711"
  //     "0383231427603";
  // Bint b1 = "3", b2 = "1928";

  // Bint b;

  // for (int i = 0; i < 100000; i++) {
  //   b1 += 1;
  //   b2 += 1;
  //   b = b1 * b2;
  // }
  // b = b1 * b2;

  // cout << b1 << endl << b2 << endl;

  // cout << b << " " << b1.bit_length() << endl;

  // Bint b3 = "1099511627776123456789";
  // pair<Bint, Bint> p = Bint::div(b3, 1000000000);
  // cout << p.first << " " << p.second << " " << b3 << endl;

  // Bint b0 =
  //     "000000077585238500233229835691539745687222933427537611140338446921599329"
  //     "563854181639226164381928796979266338729820793681714036123333769079646600"
  //     "909317387199098296728849226170703244369279252719528087576447212282120455"
  //     "268797571482789028589388060453481530946636983169469071820278066683677230"
  //     "258285779258489871121164383949525509812282559535360202183909040927935655"
  //     "246154667554657736332689570051822428740409972317367200374311517418025000"
  //     "244750063441338782593498649933548057189918817813209444543818566261537270"
  //     "164945825026980348078479972315666739582243205039432439292161591287230860"
  //     "818224740346890088322002354285178080966626692376142753101557760875126864"
  //     "991994508460081430206315243360110400087920782365366318801751222527843865"
  //     "216661603830342889155885053507741760582801868400807209962557789578344794"
  //     "184118928941770487755880556571824641816526811148918733248819867724136758"
  //     "513800381975727005271555515327063595428729629128111403012652961";

  // cout << Bint::integral_rth_root(b0, 241) << endl;

  // Bint rmax = "10000000000";
  // for (int i = 0; i < 10; i++) {
  //   cout << Bint::rand(rmax) << endl;
  // }
  // for (int i = 0; i < 10; i++) {
  //   cout << Bint::rand(-rmax) << endl;
  // }

  // Bint b1 =
  // "981273912739176243876298361982746381972683471682376182736182736",
  //      b2 = "9182739182738712638761238";

  // pair<Bint, Bint> p = Bint::div(b1, b2);

  // cout << p.first << " " << p.second << endl;

  // Bint p = "19736887067";
  // for (int i = 0; i < 10000; i++) {
  //   cout << Bint::pow(i, p - 1, p) << endl;
  // }

  //   for (Bint i =
  //            "1234567911234567911234567911234567911234567911234567911234567911234"
  //            "5679112345679112345679112345679112345679112345679112345679112345679"
  //            "1123456791123456791123456791123456791123456791123456859";
  //        ; i++) {
  //     if (Algorithm::is_probable_prime(i)) {
  //       cout << "Next prime: " << i << endl;
  //       break;
  //     }
  //   }
  // Bint i =
  //     "432243812959272959023682044933847714492379681406129049042740760037008001"
  //     "432243812959272959023682044933847714492379681406129049042740760037008001"
  //     "354678327854198237912873981273981723987128973618723681726387162387162387"
  //     "612983719283791827391827398127019238791827398172398172398712983719283793"
  //     "191209831092123123121982379128739817239812391239861298735189723891623098"
  //     "761287936891726398712638312312335467832785419823791287398127398172398712"
  //     "897361872368172638716238716238735467832712312354124381723987128973123618"
  //     "723681726387123121623871623877162387162387354678327854198237912873981273"
  //     "98712837128973618723681726387162387162387638716238716238723";
  // cout << i << endl;
  // cout << i * i << endl;

  // cout << b1 << " " << b2 << endl << b1 - b2 << endl << b2 - b1 << endl;
  // cout << Bint::gcd(696969, 69264);

  vector<vector<uint32_t>> m = {
      {1, 1, 0, 0}, {1, 1, 0, 1}, {0, 0, 1, 1}, {0, 0, 1, 0}, {0, 0, 0, 1}};

  vector<vector<int>> m2 = Algorithm::gf2_gaussian_elimination(m);

  for (auto v : m2) {
    for (auto i : v) {
      cout << i << " ";
    }
    cout << endl;
  }

  // int N = 3731;

  // vector<uint32_t> v = Algorithm::primes_less_than(N);

  // cout << v.size() << endl;
  // for (int i : v) {
  //   cout << i << " ";
  // }
  // cout << endl;

  // cout << Algorithm::square_root_modulo_prime_power(56, 101, 3) << endl;

  // Bint a = "12512312312312312312321000";
  // Bint b = "-1";
  // cout << a * b << endl;
}
