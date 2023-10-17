import java.math.BigInteger;

public class Test {
    private static BigInteger fact(BigInteger n) {
        BigInteger ans = new BigInteger("1");

        for (BigInteger i = BigInteger.TWO; i.compareTo(n) <= 0; i = i.add(BigInteger.ONE)) {
            ans = ans.multiply(i);
        }
        return ans;
    }

    private static BigInteger[] rd(BigInteger n) {
        for (int r = 0;; r++) {
            if (n.mod(BigInteger.TWO).equals(BigInteger.ONE)) return new BigInteger[] { BigInteger.valueOf(r), n };
            n = n.divide(BigInteger.TWO);
        }
    }

    public static void main(String args[]) {
        System.out.println("asdf");
        System.out.println(Test2.foo());

        final long startTime = System.currentTimeMillis();

        BigInteger a = new BigInteger("52300");
        BigInteger fac = fact(a);
        //BigInteger[] rdRes = rd(fac);

        final long endTime = System.currentTimeMillis();

        //System.out.println(rdRes[0]);
        System.out.println("Total execution time: " + (endTime - startTime) / 1000.0 + "ms");
    }
}
