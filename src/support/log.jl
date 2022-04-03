function log()

    println("\n", "="^99)

    println("Started $(StackTraces.stacktrace()[2].func)\n")

    println(now())

    println("="^99, "\n")

end
