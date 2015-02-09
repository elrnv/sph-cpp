
inc = 0.05;
count = 0;
for x = (-1+inc):inc:(1-inc)
  for y = 0+inc:inc:(1-inc)
    for z = (-1+inc):inc:(1-inc)
      printf("v %f %f %f\n", x, y, z);
      count = count + 1;
    end
  end
end

for i = 1:count
  printf("f %d\n", i);
end

