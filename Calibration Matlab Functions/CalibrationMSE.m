function arg = CalibrationMSE(u_test, v_test, u_guess, v_guess)

arg = sum(sum((u_test-u_guess).^2 + (v_test-v_guess).^2 ));

end