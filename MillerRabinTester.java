import java.math.BigInteger;
import java.util.*;
import java.math.*;

public class MillerRabinTester 
{
	public static void main(String args[])
	{
		//System.out.print(BigInteger.ZERO.compareTo(new BigInteger("9645123847672438956787934267089253")));
		//System.out.println((new BigInteger("2")).pow(511));
		//System.out.println((new BigInteger("2")).pow(512));
		//System.out.println("The random BigInteger = " + buildNBitRandom(512));
		//getBigIntFromBitString("1111");
		//randomBigInt(new BigInteger("6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048"),
		//new BigInteger("13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095"));
		//System.out.println(millerRabin(BigInteger.valueOf(65537), BigInteger.valueOf(10)));
		//boolean mayBePrime = false;
		/*while(!mayBePrime)
		{
			mayBePrime = millerRabin(randomBigInt(new BigInteger("6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048"),
					new BigInteger("13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095")), BigInteger.valueOf(10));
		}*/
		/*boolean isStrongPrime = false;
		while(!isStrongPrime)
		{
			isStrongPrime = true;
		}*/
		//dumbFactoring(BigInteger.valueOf(65536));
		//System.out.println(new BigInteger("18446744073709551617").divide(BigInteger.valueOf(274177)));
		//dumbFactoring(new BigInteger("18446744073709551617"));
		
		pollardRho(new BigInteger("18446744073709551617"));
		//System.out.println((new BigInteger("2")).pow(128));
		//System.out.println((new BigInteger("340282366920938463463374607431768211457").divide(new BigInteger("59649589127497217"))));
	}

	public static BigInteger[] hasLargePrimeFactors(BigInteger p)
	{
		BigInteger[] arr = {BigInteger.valueOf(0),BigInteger.valueOf(0)};

		return arr;
	}

	public static ArrayList<BigInteger> dumbFactoring(BigInteger x)
	{
		ArrayList<BigInteger> list = new ArrayList<BigInteger>();
		BigInteger n = BigInteger.valueOf(2);
		BigInteger temp = x;
		while(n.compareTo(x.sqrt()) == -1 && !temp.equals(BigInteger.ONE))
		{
			if(temp.mod(n).equals(BigInteger.ZERO))
			{
				temp = temp.divide(n);
				list.add(n);
			}
			else 
				n = n.add(BigInteger.ONE);
			System.out.println(n);
		}
		System.out.println(list.toString());
		return list;
	}

	public static BigInteger randomBigInt(BigInteger minLimit, BigInteger maxLimit)
	{
		BigInteger bigInteger = maxLimit.subtract(minLimit);
		Random randNum = new Random();
		int len = maxLimit.bitLength();
		BigInteger res = new BigInteger(len, randNum);
		if (res.compareTo(minLimit) < 0)
			res = res.add(minLimit);
		if (res.compareTo(bigInteger) >= 0)
			res = res.mod(bigInteger).add(minLimit);
		//System.out.println("The random BigInteger = "+ res);
		return res;
	}

	/**
	 * @param l length of the integer to be generated, in bits
	 * 
	 * @return res, the BigInteger in decimal notation
	 */
	public static BigInteger buildNBitRandom(int l)
	{
		//hold a temporary int that we randomly generate
		int temp = 0; 
		//make sure there is a leading 1 and l - 1 bits to the right
		String num = "1";
		/*
		 * This method works by generating a random int, then converting it to binary,
		 * then concatenating it's bitstring representation to "num" until num is of length
		 * at least l
		 */
		while(num.length() < l)
		{
			temp = (int) (Math.random() * Integer.MAX_VALUE);
			num = num + Integer.toBinaryString(temp);
			//System.out.println(Integer.toBinaryString(temp));
		}
		//System.out.println(num.substring(0,l));
		//System.out.println(num.substring(0, l).length());
		//once num has length >= l then we can turn it back from binary to BigInteger
		BigInteger res = getBigIntFromBitString(num.substring(0,l));
		//System.out.println(res);
		return res;
	}

	/**
	 * 
	 * @param s the string we want to convert into a BigInteger
	 * 
	 * @return the string s as a BigInteger representation
	 */
	public static BigInteger getBigIntFromBitString(String s)
	{
		BigInteger res = new BigInteger("0");
		BigInteger two = new BigInteger("2");
		//for every column take the bit and check if its a 1 and if it is add 2^n based on the row
		for(int i = 0; i < s.length(); i++)
		{
			if(Character.getNumericValue(s.charAt(i)) == 1)
			{	
				res = res.add(two.pow(s.length() - i - 1));
			}
		}
		//System.out.println(res);
		return res;
	}

	/**
	 * 
	 * @param n integer we want to test
	 * @param k number of rounds to do 
	 * @return true if n may be prime, false if it is not prime
	 */
	public static boolean millerRabin(BigInteger n, BigInteger k)
	{
		BigInteger s = BigInteger.ZERO;
		BigInteger d = BigInteger.ZERO;
		BigInteger two = BigInteger.valueOf(2);
		//n - 1 which we must factor by 2 in order to find d and s
		BigInteger temp = n.subtract(BigInteger.ONE);
		BigInteger a = BigInteger.ZERO;
		BigInteger x = BigInteger.ZERO;
		BigInteger y = BigInteger.ZERO;
		//keep dividing temp by 2 until you cannot divide exactly, while increasing 
		//s every time that you divide by 2 and get an integer
		//then change your d to whatever remains
		while(d.equals(BigInteger.ZERO))
		{
			if(temp.mod(two).equals(BigInteger.ZERO))
			{
				s = s.add(BigInteger.ONE);
				temp = temp.divide(two);
			}
			else 
				d = temp;
		}
		//System.out.println(n.subtract(BigInteger.ONE) + " = 2^" + s + " * " + d);
		//"repeat k times"
		for(BigInteger i = BigInteger.ZERO; i.compareTo(k) < 0; i = i.add(BigInteger.ONE))
		{
			a = randomBigInt(BigInteger.valueOf(2), n.subtract(BigInteger.valueOf(2)));
			x = a.modPow(d, n);
			//"repeat s times"
			for(BigInteger j = BigInteger.ZERO; j.compareTo(s) < 0; j = j.add(BigInteger.ONE))
			{
				y = x.modPow(BigInteger.valueOf(2), n);
				if(y.equals(BigInteger.ONE) && !(x.equals(BigInteger.ONE)) && !(x.equals(n.subtract(BigInteger.ONE))))
					return false;
				x = y;
			}
			if(!(y.equals(BigInteger.ONE)))
				return false;
		}
		System.out.println(n);
		return true;
	}

	//finds a^-1 mod n or ax congruent to 1 mod n
	/*boolean isGcdOne(BigInteger a, BigInteger n)
	{
		//These arrays store the values of each step of the euclidean algorithm as a 4-tuple
		//this allows us to recover them when we want to find the value X that multiplies to 1 mod n
		/*
		 * For example let's say we had gcd(6,17) then origNums[0] = 17 nums[0] = 6 factors[0] = 2 and remainders[0] = 5
		 * and follow this until our remainders[i] is 1
		 
		ArrayList<BigInteger> nums = new ArrayList<BigInteger>();
		ArrayList<BigInteger> origNums = new ArrayList<BigInteger>();
		ArrayList<BigInteger> factors = new ArrayList<BigInteger>();
		ArrayList<BigInteger> remainders = new ArrayList<BigInteger>();
		//the value that we start at 
		BigInteger valueToBreakUp = n;
		//see how many times this value fits into "valueToBreakUp"
		BigInteger rightVal = a;
		//how many times rightVal fits into valueToBreakUp
		BigInteger factor = n.divide(a);
		//remainder of the division
		BigInteger remainder = n.mod(a);
		//used to count how many steps it took to find the gcd of the 2 numbers
		BigInteger i = BigInteger.ONE;
		//later becomes equal to the amount of steps we took and its used in order to find the inverse
		BigInteger lastStep = BigInteger.ZERO;
		//store our 4 tuple for a given step where origNums[i] = nums[i] * factors[i] + remainders[i]
		origNums.add(valueToBreakUp);
		nums.add(rightVal);
		factors.add(factor);
		remainders.add(remainder);

		//print the tuple to see the whole algorithm
		//printf("%lf = %lf * %lf + %lf\n", valueToBreakUp, rightVal, factor, remainder);
		
		//now iterate backwards to find all the equivalences such that origNums[i] - nums[i] * factors[i] = remainders[i]
		while (remainder.compareTo(BigInteger.ONE) > 0 && i.compareTo(n) == -1)
		{
			valueToBreakUp = rightVal;
			rightVal = remainder;
			factor = valueToBreakUp.divide(rightVal).subtract(BigInteger.ONE);
			remainder = valueToBreakUp.mod(rightVal);
			//printf("%lf = %lf * %lf + %lf\n", valueToBreakUp, rightVal, factor, remainder);
			origNums.add(valueToBreakUp);
			nums.add(rightVal);
			factors.add(factor);
			remainders.add(remainder);
			//printf("%lf + %lf * %lf = %lf\n", valueToBreakUp, rightVal, -factor, remainder);
			i.add(BigInteger.ONE);
		}
		lastStep = i;
		i = BigInteger.ZERO;

		//these variables are used to store all the values of a tuple and do the calculations to transform them over 
		//to the next step while mantaining the same structure of k1 * some factor + k2 * some factor = 1 
		if (lastStep > 1)
		{

			double frontFactor = 1.0;
			double backFactor = 0.0;
			if (factors)
				backFactor = -factors[lastStep - 1];
			double frontNum = 0.0;
			if (origNums)
				frontNum = origNums[lastStep - 1];
			double backNum = 0.0;
			if (nums)
				backNum = nums[lastStep - 1];
			double congruenceRemainder = 1.0;
			if (remainders)
				congruenceRemainder = remainders[lastStep - 1];
			//printf("(%lf) %lf + %lf (%lf) = %lf\n", frontFactor, frontNum, backNum, backFactor, congruenceRemainder);
			//do this until we have reached the last of our tuples and we have n * some factors + a * some factor = 1 mod n
			while (lastStep - 1 > 0)
			{
				//the process will basically alternate based on which of the 2 numbers is smaller, first the right will be smaller so go into this step
				//then the left will be smaller and so we will use the "else" part of this loop, this method works only if the gcd of the numbers is 1 
				if ((lastStep - 1) % 2 == 0)
				{
					frontFactor += frontFactor * -backFactor;
					backNum = origNums[lastStep - 2];
					backFactor = -factors[lastStep - 2];
					//printf("(%lf) %lf + %lf (%lf) = %lf\n", frontFactor, frontNum, backNum, backFactor, congruenceRemainder);
				}
				else
				{
					if (factors)
						backFactor = frontFactor * -factors[lastStep - 2] + backFactor;
					if (nums)
						backNum = nums[lastStep - 2];
					if (origNums)
						frontNum = origNums[lastStep - 2];
					//printf("(%lf) %lf + %lf (%lf) = %lf\n", frontFactor, frontNum, backNum, backFactor, congruenceRemainder);
				}
				lastStep--;
			}

			//printf("%lf\n", -backFactor);
			//since we are working with powers of the inverse when using 
			//Shanks' algorithm, we rather have a positive number for the inverse
			if (-backFactor < 0)
				return -backFactor + n;
			return -backFactor;
		}
		return factor;
	}*/
	
    public static BigInteger gcd(BigInteger a, BigInteger b)
    {
        // if b=0, a is the GCD
        if (b.equals(BigInteger.ZERO))
            return a;
 
        // call the gcd() method recursively by
        // replacing a with b and b with
        // modulus(a,b) as long as b != 0
        else
            return gcd(b, a.mod(b));
    }
	
	public static ArrayList<BigInteger> pollardRho(BigInteger n)
	{
		ArrayList<BigInteger> list = new ArrayList<BigInteger>();
		BigInteger x = BigInteger.valueOf(2);
		BigInteger y = BigInteger.valueOf(2);
		BigInteger d = BigInteger.ONE;
		while(d.equals(BigInteger.ONE))
		{
			x = x.pow(2).add(BigInteger.ONE).mod(n);
			y = ((y.pow(2).add(BigInteger.ONE).mod(n)).pow(2)).add(BigInteger.ONE).mod(n);
			d = gcd(x.subtract(y).abs(), n);
			System.out.println("x = " + x + " y = " + y);
		}
		if(d.equals(n))
			return list;
		list.add(d);
		System.out.println(d);
		return list;
	}
}
